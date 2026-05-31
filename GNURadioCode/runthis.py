#!/usr/bin/env python3
"""
analyze_shift.py  (peak-tracking + before/after visualization)

Measures the average before->after frequency SHIFT and peak BROADENING of RFI
peaks between a stationary ("before") and rotating ("after") recording, and
renders a presentable before/after figure.

Method (generalises Figure 3 of the DAWC paper):
  1. Load both spectra.
  2. Detect strong RFI peaks ONCE, in the stationary spectrum.
  3. For each, locate the same peak in a small window of the rotating spectrum,
     refined with quadratic (sub-bin) interpolation -> SHIFT.
  4. Compare second-moment widths before vs after -> BROADENING.
  5. Average across peaks; judge against a stationary-vs-stationary null floor.

    python3 analyze_shift.py stationary.dat rotating.dat
    python3 analyze_shift.py stationary.dat rotating.dat --plot shift.png
"""

import argparse
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d

import matplotlib
matplotlib.use("Agg")                      # render to file without a display
import matplotlib.pyplot as plt
from matplotlib import gridspec


def load_spectrum(path, fft_size, mode='last'):
    data = np.fromfile(path, dtype=np.float32)
    n = data.size // fft_size
    if n == 0:
        raise SystemExit(f"{path}: file too short / wrong fft-size.")
    frames = data[:n * fft_size].reshape(n, fft_size)
    return frames[-1] if mode == 'last' else frames.mean(axis=0)


def snr_db_to_height(snr_db, sigma):
    return max(snr_db, 4.0 * sigma)


def detect_reference_peaks(spec, prominence, snr_db, top_n):
    """Detect only genuine RFI peaks well above the noise floor."""
    med = np.median(spec)
    mad = np.median(np.abs(spec - med)) + 1e-9
    sigma = 1.4826 * mad
    height = med + snr_db_to_height(snr_db, sigma)
    peaks, _ = find_peaks(spec, prominence=prominence, height=height)
    if peaks.size == 0:
        peaks, _ = find_peaks(spec, prominence=prominence, height=med + 3 * sigma)
    if peaks.size == 0:
        raise SystemExit("No peaks rise above the noise floor; "
                         "lower --snr or check that real RFI is present.")
    if peaks.size > top_n:
        order = np.argsort(spec[peaks])[::-1][:top_n]
        peaks = np.sort(peaks[order])
    return peaks


def refine_peak(spec, k):
    """3-point parabolic interpolation -> fractional bin index (sub-bin)."""
    if k <= 0 or k >= len(spec) - 1:
        return float(k)
    a, b, c = spec[k - 1], spec[k], spec[k + 1]
    denom = a - 2.0 * b + c
    if denom == 0:
        return float(k)
    return k + 0.5 * (a - c) / denom


def rms_width_var(spec_db, center_bin, half_win, floor_db):
    """Intensity-weighted variance (bins^2) of a peak, on linear power with the
    noise floor subtracted."""
    lo = max(0, center_bin - half_win)
    hi = min(len(spec_db), center_bin + half_win + 1)
    idx = np.arange(lo, hi)
    p = 10.0 ** (spec_db[lo:hi] / 10.0) - 10.0 ** (floor_db / 10.0)
    p = np.clip(p, 0.0, None)
    tot = p.sum()
    if tot <= 0:
        return np.nan
    centroid = (p * idx).sum() / tot
    return float((p * (idx - centroid) ** 2).sum() / tot)


def make_plot(stat, rot, records, bw, center_freq, fft_size, summary,
              outpath, n_panels):
    """Overview of the whole spectrum + zoomed before/after panels."""
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    c_stat, c_rot = "#2c6fbb", "#e8743b"
    n = len(stat)
    f_mhz = (center_freq + (np.arange(n) - fft_size / 2) * bw) / 1e6

    recs = sorted(records, key=lambda r: r['power'], reverse=True)[:n_panels]
    recs = sorted(recs, key=lambda r: r['sp'])
    ncols = max(len(recs), 1)

    fig = plt.figure(figsize=(3.3 * max(ncols, 3), 7.4))
    gs = gridspec.GridSpec(2, ncols, height_ratios=[1.05, 1.0],
                           hspace=0.42, wspace=0.30)

    # ---- overview ----
    ax0 = fig.add_subplot(gs[0, :])
    ax0.plot(f_mhz, stat, color=c_stat, lw=0.8, alpha=0.9, label='Stationary (before)')
    ax0.plot(f_mhz, rot, color=c_rot, lw=0.8, alpha=0.7, label='Rotating (after)')
    for r in records:
        b = int(round(r['sp']))
        ax0.plot(f_mhz[b], stat[b], 'v', color=c_stat, ms=6, mec='white', mew=0.6)
    ax0.set_title('Full spectrum — tracked RFI peaks (▼)', fontweight='bold')
    ax0.set_xlabel('Frequency (MHz)')
    ax0.set_ylabel('Power (dB)')
    ax0.legend(loc='upper right', framealpha=0.9)
    ax0.margins(x=0.01)

    # ---- zoom panels ----
    for i, r in enumerate(recs):
        ax = fig.add_subplot(gs[1, i])
        lo, hi = r['lo'], r['hi']
        x = (np.arange(lo, hi) - r['sp_ref']) * bw           # Hz rel. to before-peak
        ax.plot(x, stat[lo:hi], '-o', color=c_stat, ms=3, lw=1.3, label='before')
        ax.plot(x, rot[lo:hi], '-o', color=c_rot, ms=3, lw=1.3, label='after')
        shift = (r['rp_ref'] - r['sp_ref']) * bw
        ax.axvline(0, color=c_stat, ls='--', lw=1)
        ax.axvline(shift, color=c_rot, ls='--', lw=1)
        ax.set_title(f"{f_mhz[int(round(r['sp']))]:.3f} MHz\n"
                     f"$\\Delta f$ = {shift:+.0f} Hz", fontsize=9)
        ax.set_xlabel('$\\Delta f$ from before-peak (Hz)', fontsize=8)
        if i == 0:
            ax.set_ylabel('Power (dB)')
        ax.legend(fontsize=7, loc='upper right')
        ax.tick_params(labelsize=8)

    fig.suptitle(summary, fontsize=12.5, fontweight='bold', y=0.995)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="RFI peak shift + broadening, with plot.")
    p.add_argument('stationary')
    p.add_argument('rotating')
    p.add_argument('--fft-size', type=int, default=8192)
    p.add_argument('--samp-rate', type=float, default=2.4e6)
    p.add_argument('--center-freq', type=float, default=1.42041e9)
    p.add_argument('--prominence', type=float, default=6.0)
    p.add_argument('--snr', type=float, default=8.0)
    p.add_argument('--top-n', type=int, default=40)
    p.add_argument('--window-hz', type=float, default=1500.0)
    p.add_argument('--width-hz', type=float, default=2500.0)
    p.add_argument('--average', choices=['last', 'mean'], default='last')
    p.add_argument('--smooth', type=int, default=3)
    p.add_argument('--plot', default='shift_plot.png',
                   help='output figure path (PNG)')
    p.add_argument('--plot-peaks', type=int, default=4,
                   help='how many of the strongest peaks to show as zoom panels')
    args = p.parse_args()

    bw = args.samp_rate / args.fft_size
    W = max(2, int(round(args.window_hz / bw)))
    Wd = max(3, int(round(args.width_hz / bw)))

    stat = load_spectrum(args.stationary, args.fft_size, args.average)
    rot = load_spectrum(args.rotating, args.fft_size, args.average)
    if args.smooth and args.smooth > 1:
        stat = uniform_filter1d(stat, size=args.smooth)
        rot = uniform_filter1d(rot, size=args.smooth)
    n = len(stat)
    floor_db = float(np.median(stat))

    ref_peaks = detect_reference_peaks(stat, args.prominence, args.snr, args.top_n)

    shifts_hz, weights, var_stat, var_rot, records = [], [], [], [], []
    for sp in ref_peaks:
        lo, hi = max(0, sp - W), min(n, sp + W + 1)
        if hi - lo < 3:
            continue
        rp = lo + int(np.argmax(rot[lo:hi]))
        sp_ref = refine_peak(stat, sp)
        rp_ref = refine_peak(rot, rp)
        shifts_hz.append((rp_ref - sp_ref) * bw)
        weights.append(stat[sp])
        vs = rms_width_var(stat, sp, Wd, floor_db)
        vr = rms_width_var(rot, sp, Wd, floor_db)
        if np.isfinite(vs) and np.isfinite(vr):
            var_stat.append(vs)
            var_rot.append(vr)
        records.append({'sp': sp, 'sp_ref': sp_ref, 'rp_ref': rp_ref,
                        'lo': lo, 'hi': hi, 'power': stat[sp]})

    shifts = np.array(shifts_hz)
    mag = np.abs(shifts)
    w = np.array(weights)
    w = w - w.min() + 1e-9

    mean_vs = float(np.nanmean(var_stat)) if var_stat else float('nan')
    mean_vr = float(np.nanmean(var_rot)) if var_rot else float('nan')
    doppler_var = mean_vr - mean_vs
    doppler_rms_hz = (doppler_var ** 0.5 * bw) if doppler_var > 0 else 0.0

    print(f"Reference peaks tracked:     {len(shifts)}  (strong RFI only)")
    print(f"Frequency resolution:        {bw:.1f} Hz/bin")
    print("--- shift metric (peak displacement) ---")
    print(f"Mean |shift|:                {mag.mean():8.1f} Hz   <-- headline")
    print(f"Power-weighted mean |shift|: {np.average(mag, weights=w):8.1f} Hz")
    print(f"Median |shift|:              {np.median(mag):8.1f} Hz")
    print(f"Mean signed shift:           {shifts.mean():+8.1f} Hz   (drift/bias check)")
    print("--- broadening metric (peak width growth) ---")
    print(f"Mean width, stationary:      {mean_vs ** 0.5 * bw:8.1f} Hz (RMS)")
    print(f"Mean width, rotating:        {mean_vr ** 0.5 * bw:8.1f} Hz (RMS)")
    print(f"Doppler spread (added RMS):  {doppler_rms_hz:8.1f} Hz   <-- relative indicator")

    summary = (f"RFI peaks: stationary vs rotating    |    "
               f"mean |shift| {mag.mean():.0f} Hz   •   "
               f"Doppler spread {doppler_rms_hz:.0f} Hz   •   {len(shifts)} peaks")
    make_plot(stat, rot, records, bw, args.center_freq, args.fft_size,
              summary, args.plot, args.plot_peaks)
    print(f"\nFigure written to {args.plot}")


if __name__ == '__main__':
    main()