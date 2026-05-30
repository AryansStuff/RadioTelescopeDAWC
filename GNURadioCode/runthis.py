#!/usr/bin/env python3
"""
analyze_shift.py  (average-over-peaks version)

Generalises the paper's measurement (Figure 3). Instead of using only the single
highest peak, it matches EVERY peak in the stationary ("before") spectrum to its
nearest peak in the rotating ("after") spectrum, measures each before->after
frequency distance with quadratic interpolation, and averages them so the result
represents the spectrum as a whole.

    python3 analyze_shift.py stationary.dat rotating.dat

The dominant-peak number from the original method is still printed for continuity.
"""

import argparse
import numpy as np
from scipy.signal import find_peaks


def load_spectrum(path, fft_size, mode='last'):
    """Return one 8192-bin spectrum: the IIR-converged last frame, or the
    mean of every saved frame."""
    data = np.fromfile(path, dtype=np.float32)
    n = data.size // fft_size
    if n == 0:
        raise SystemExit(f"{path}: file too short / wrong fft-size.")
    frames = data[:n * fft_size].reshape(n, fft_size)
    return frames[-1] if mode == 'last' else frames.mean(axis=0)


def refine_peak(spec, k):
    """3-point parabolic interpolation -> fractional bin index."""
    if k <= 0 or k >= len(spec) - 1:
        return float(k)
    a, b, c = spec[k - 1], spec[k], spec[k + 1]
    denom = a - 2.0 * b + c
    if denom == 0:
        return float(k)
    return k + 0.5 * (a - c) / denom


def main():
    p = argparse.ArgumentParser(description="Average RFI peak shift across all peaks.")
    p.add_argument('stationary', help='stationary recording (e.g. stationary.dat)')
    p.add_argument('rotating', help='rotating recording (e.g. rotating.dat)')
    p.add_argument('--fft-size', type=int, default=8192)
    p.add_argument('--samp-rate', type=float, default=2.4e6)
    p.add_argument('--center-freq', type=float, default=1.42041e9)
    p.add_argument('--prominence', type=float, default=3.0,
                   help='min peak prominence (dB) for a bin to count as a peak')
    p.add_argument('--max-shift', type=float, default=2000.0,
                   help='largest plausible shift (Hz) when matching a before peak '
                        'to an after peak; prevents matching to a different source')
    p.add_argument('--average', choices=['last', 'mean'], default='last')
    args = p.parse_args()

    bw = args.samp_rate / args.fft_size            # Hz per bin
    max_bins = max(1, int(round(args.max_shift / bw)))

    stat = load_spectrum(args.stationary, args.fft_size, args.average)
    rot = load_spectrum(args.rotating, args.fft_size, args.average)

    stat_peaks, _ = find_peaks(stat, prominence=args.prominence)
    rot_peaks, _ = find_peaks(rot, prominence=args.prominence)
    if stat_peaks.size == 0 or rot_peaks.size == 0:
        raise SystemExit("Not enough peaks found; lower --prominence.")

    shifts_hz, weights = [], []
    for sp in stat_peaks:
        d = np.abs(rot_peaks - sp)                 # nearest rotating peak
        j = int(np.argmin(d))
        if d[j] > max_bins:
            continue                               # no plausible match -> skip
        sp_ref = refine_peak(stat, sp)
        rp_ref = refine_peak(rot, rot_peaks[j])
        shifts_hz.append((rp_ref - sp_ref) * bw)
        weights.append(stat[sp])                   # peak power as a weight

    if not shifts_hz:
        raise SystemExit("No peaks matched within --max-shift; widen it.")

    shifts = np.array(shifts_hz)
    mag = np.abs(shifts)
    w = np.array(weights)
    w = w - w.min() + 1e-9                          # keep weights positive

    print(f"Matched peaks:               {len(shifts)} of {len(stat_peaks)} stationary peaks")
    print(f"Frequency resolution:        {bw:.1f} Hz/bin")
    print(f"Mean |shift|:                {mag.mean():8.1f} Hz   <-- headline (scales with RPM)")
    print(f"Median |shift|:              {np.median(mag):8.1f} Hz")
    print(f"Power-weighted mean |shift|: {np.average(mag, weights=w):8.1f} Hz")
    print(f"Std of |shift|:              {mag.std():8.1f} Hz")
    print(f"Mean signed shift:           {shifts.mean():+8.1f} Hz   (net drift / bias check)")

    # Original single-dominant-peak result, kept for continuity with the paper
    kr = int(np.argmax(rot))
    kr_ref = refine_peak(rot, kr)
    near = stat_peaks[int(np.argmin(np.abs(stat_peaks - kr)))]
    ks_ref = refine_peak(stat, near)
    print(f"Dominant-peak shift:         {(kr_ref - ks_ref) * bw:+8.1f} Hz")


if __name__ == '__main__':
    main()