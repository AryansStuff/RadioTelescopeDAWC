#!/usr/bin/env python3
"""
This needs RAW IQ, not averaged spectra. Record it with this same script:

    python3 doppler_waterfall.py record rotating.bin --center 1420.4e6 --rate 2.4e6 --seconds 30 --gain 40 --warmup 60

You can also record and analyze in one shot:

    python3 doppler_waterfall.py record rotating.bin --center 1420.4e6 --rate 2.4e6 --seconds 30 --analyze --rpm 100 --radius 1.0

ANALYSIS
--------
    python3 doppler_waterfall.py analyze rotating.bin --rate 2.4e6 --center 1420.4e6 --rpm 100 --radius 1.0

  # if you recorded 8-bit IQ from rtl_sdr instead: --format u8
  # pick a specific RFI peak:                       --peak-offset 150000
"""

import argparse
import sys

import numpy as np

try:
    from scipy import signal
except ImportError:
    sys.exit("Missing dependency. Run:  pip install numpy scipy matplotlib")

import matplotlib
matplotlib.use("Agg")  # headless: save PNG without a display
import matplotlib.pyplot as plt

C = 299_792_458.0  # speed of light, m/s


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_iq(path, fmt, rate, max_seconds):
    """Load complex IQ samples, capped at max_seconds to bound memory."""
    max_samp = int(max_seconds * rate) if max_seconds else None
    if fmt == "cf32":  # GNU Radio gr_complex = interleaved float32
        data = np.fromfile(path, dtype=np.complex64, count=(max_samp or -1))
        x = data
    elif fmt == "u8":  # rtl_sdr raw: interleaved unsigned bytes, 127.5-centered
        raw = np.fromfile(path, dtype=np.uint8,
                          count=(2 * max_samp if max_samp else -1))
        raw = raw[: 2 * (len(raw) // 2)]
        x = ((raw[0::2].astype(np.float32) - 127.5) +
             1j * (raw[1::2].astype(np.float32) - 127.5)) / 127.5
        x = x.astype(np.complex64)
    elif fmt == "i16":  # interleaved signed 16-bit
        raw = np.fromfile(path, dtype=np.int16,
                          count=(2 * max_samp if max_samp else -1))
        raw = raw[: 2 * (len(raw) // 2)]
        x = (raw[0::2].astype(np.float32) +
             1j * raw[1::2].astype(np.float32)) / 32768.0
        x = x.astype(np.complex64)
    else:
        sys.exit(f"Unknown --format {fmt!r} (use cf32, u8, or i16)")
    if x.size == 0:
        sys.exit("No samples loaded - check the file path and --format.")
    return x


# ---------------------------------------------------------------------------
# Find the RFI peak in the band
# ---------------------------------------------------------------------------
def find_peak_offset(x, rate, nfft=8192, dc_guard_hz=3000.0):
    """Average power spectrum; return the offset (Hz) of the strongest peak,
    ignoring a small guard band around DC (the RTL-SDR's center spike)."""
    f, pxx = signal.welch(x, fs=rate, nperseg=min(nfft, x.size),
                          return_onesided=False, scaling="density")
    f = np.fft.fftshift(f)
    pxx = np.fft.fftshift(pxx)
    mask = np.abs(f) > dc_guard_hz
    idx = np.argmax(pxx[mask])
    return float(f[mask][idx])


# ---------------------------------------------------------------------------
# Channelize: shift the chosen peak to baseband and decimate for fine resolution
# ---------------------------------------------------------------------------
def _factorize(q, cap=10):
    """Break a large decimation factor into stages each <= cap."""
    factors = []
    for p in (10, 9, 8, 7, 6, 5, 4, 3, 2):
        while q % p == 0 and p <= cap:
            factors.append(p)
            q //= p
    if q > 1:
        factors.append(q)  # leftover (decimate handles moderate factors)
    return factors


def channelize(x, rate, peak_offset, target_rate):
    """Mix the peak to 0 Hz and decimate to ~target_rate so the spectrogram is
    small and high-resolution. Returns (baseband_signal, new_rate)."""
    n = np.arange(x.size, dtype=np.float64)
    x0 = x * np.exp(-2j * np.pi * peak_offset * n / rate).astype(np.complex64)
    decim = max(1, int(round(rate / target_rate)))
    if decim > 1:
        for q in _factorize(decim):
            if q > 1:
                x0 = signal.decimate(x0, q, ftype="fir", zero_phase=True)
        rate = rate / decim
    return x0.astype(np.complex64), rate


# ---------------------------------------------------------------------------
# Waterfall + per-slice peak tracking
# ---------------------------------------------------------------------------
def _parabolic(ym1, y0, yp1):
    """Sub-bin peak offset (in bins) from a 3-point parabola fit."""
    denom = (ym1 - 2 * y0 + yp1)
    return 0.0 if denom == 0 else 0.5 * (ym1 - yp1) / denom


def spectrogram_and_track(xc, rate, slice_sec, search_hz):
    """Build the spectrogram and track the peak (near 0 Hz) in each time slice."""
    nper = max(16, int(slice_sec * rate))
    noverlap = nper // 2
    f, t, Sxx = signal.spectrogram(xc, fs=rate, nperseg=nper, noverlap=noverlap,
                                   return_onesided=False, mode="psd")
    f = np.fft.fftshift(f)
    Sxx = np.fft.fftshift(Sxx, axes=0)

    in_win = np.abs(f) <= search_hz
    win_idx = np.where(in_win)[0]
    df = f[1] - f[0]

    track = np.full(Sxx.shape[1], np.nan)
    for j in range(Sxx.shape[1]):
        col = Sxx[win_idx, j]
        k = int(np.argmax(col))
        if 0 < k < len(col) - 1:
            delta = _parabolic(col[k - 1], col[k], col[k + 1])
        else:
            delta = 0.0
        track[j] = f[win_idx[0]] + (k + delta) * df
    return f, t, Sxx, track


# ---------------------------------------------------------------------------
# Fit a sine at the known rotation frequency (linear least squares)
# ---------------------------------------------------------------------------
def fit_sine(t, y, f_rot, detrend_deg=2):
    """Fit y(t) = offset + A*cos(2*pi*f_rot*t) + B*sin(2*pi*f_rot*t).
    Optionally remove a slow polynomial trend first (kills oscillator drift)."""
    good = np.isfinite(y)
    t, y = t[good], y[good]

    trend = np.zeros_like(y)
    if detrend_deg and len(y) > detrend_deg + 1:
        coef = np.polyfit(t, y, detrend_deg)
        trend = np.polyval(coef, t)
    yd = y - trend

    w = 2 * np.pi * f_rot
    M = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    (c0, A, B), *_ = np.linalg.lstsq(M, yd, rcond=None)

    fit = M @ np.array([c0, A, B])
    amplitude = float(np.hypot(A, B))
    phase = float(np.arctan2(-B, A))
    resid = yd - fit
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((yd - yd.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {
        "amplitude_hz": amplitude,
        "phase_rad": phase,
        "r2": r2,
        "resid_rms_hz": float(np.sqrt(np.mean(resid**2))),
        "t": t, "y": y, "trend": trend, "fit_detrended": fit,
    }


def track_periodogram(t, y, detrend_deg=2):
    """FFT of the (detrended) frequency-vs-time track. A real Doppler effect
    shows power concentrated at the rotation frequency."""
    good = np.isfinite(y)
    t, y = t[good], y[good]
    if detrend_deg and len(y) > detrend_deg + 1:
        y = y - np.polyval(np.polyfit(t, y, detrend_deg), t)
    dt = np.median(np.diff(t))
    Y = np.abs(np.fft.rfft(y * np.hanning(len(y))))
    fr = np.fft.rfftfreq(len(y), dt)
    return fr, Y


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def make_plot(out, f, t, Sxx, track, fitres, fr, Y, f_rot, center, peak_offset):
    abs_f = (center + peak_offset + f) / 1e6  # MHz on the waterfall axis
    fig, ax = plt.subplots(3, 1, figsize=(10, 12))

    # (1) Waterfall with tracked peak
    Sdb = 10 * np.log10(Sxx + 1e-20)
    im = ax[0].pcolormesh(t, abs_f, Sdb, shading="auto", cmap="viridis")
    ax[0].plot(t, (center + peak_offset + track) / 1e6, "r-", lw=1.0,
               label="tracked peak")
    ax[0].set_xlabel("time (s)"); ax[0].set_ylabel("frequency (MHz)")
    ax[0].set_title("Waterfall (spectrogram) with tracked RFI peak")
    ax[0].legend(loc="upper right"); fig.colorbar(im, ax=ax[0], label="dB")

    # (2) Frequency track + fitted sine
    ax[1].plot(fitres["t"], fitres["y"] - fitres["trend"], ".", ms=3,
               alpha=0.6, label="peak shift (drift removed)")
    ax[1].plot(fitres["t"], fitres["fit_detrended"], "r-", lw=2,
               label=f"sine fit @ {f_rot:.3f} Hz")
    ax[1].set_xlabel("time (s)"); ax[1].set_ylabel("frequency shift (Hz)")
    ax[1].set_title(f"Peak shift vs time  -  fitted amplitude "
                    f"{fitres['amplitude_hz']:.1f} Hz, R\u00b2={fitres['r2']:.2f}")
    ax[1].legend(loc="upper right"); ax[1].grid(alpha=0.3)

    # (3) Periodogram of the track (spike should sit at the rotation frequency)
    ax[2].plot(fr, Y, "b-")
    ax[2].axvline(f_rot, color="r", ls="--", label=f"rotation freq {f_rot:.3f} Hz")
    ax[2].set_xlim(0, max(5 * f_rot, fr[1] * 10))
    ax[2].set_xlabel("frequency (Hz)"); ax[2].set_ylabel("amplitude")
    ax[2].set_title("Periodogram of the frequency track")
    ax[2].legend(loc="upper right"); ax[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"Saved plot: {out}")


# ---------------------------------------------------------------------------
# Recording: stream raw IQ from the RTL-SDR straight to disk (memory stays flat)
# ---------------------------------------------------------------------------
def record_iq(path, rate, center, gain, seconds, warmup=0.0,
              block=256 * 1024):
    """Capture raw IQ to `path` as complex64. Writes block-by-block, so RAM use
    stays tiny no matter how long you record - this is what makes it Pi-safe."""
    try:
        from rtlsdr import RtlSdr
    except ImportError:
        sys.exit("Recording needs pyrtlsdr. Run:  pip install pyrtlsdr\n"
                 "(If the import crashes, you likely need: pip install \"pyrtlsdr==0.2.91\")")

    sdr = RtlSdr()
    sdr.sample_rate = rate
    sdr.center_freq = center
    try:
        sdr.gain = "auto" if str(gain).lower() == "auto" else float(gain)
    except Exception:
        sdr.gain = "auto"

    print(f"RTL-SDR: {rate/1e6:.3f} MS/s, center {center/1e6:.4f} MHz, gain {gain}")

    # Warm-up: read and DISCARD for a while so the oscillator settles before
    # recording (drift is worst right after power-on). Strongly recommended.
    if warmup > 0:
        print(f"Warming up {warmup:.0f} s (discarding samples to let drift settle)...")
        end = int(warmup * rate)
        got = 0
        while got < end:
            sdr.read_samples(block)
            got += block

    total = int(seconds * rate)
    written = 0
    bytes_per = 8  # complex64
    print(f"Recording {seconds:.0f} s -> {path} "
          f"(~{total*bytes_per/1e6:.0f} MB on disk)...")
    try:
        with open(path, "wb") as f:
            while written < total:
                samples = sdr.read_samples(block)
                samples.astype(np.complex64).tofile(f)
                written += len(samples)
                pct = 100.0 * written / total
                print(f"\r  {written/rate:6.1f}/{seconds:.0f} s  ({pct:4.0f}%)",
                      end="", flush=True)
        print()
    finally:
        sdr.close()
    print(f"Done. Recorded {written/rate:.1f} s = {written:,} samples.")
    return written / rate


# ---------------------------------------------------------------------------
def run_analysis(args):
    if args.rot_freq:
        f_rot = args.rot_freq
    elif args.rpm:
        f_rot = args.rpm / 60.0
    else:
        sys.exit("Provide --rpm or --rot-freq (the rotation rate to fit against).")

    rot_period = 1.0 / f_rot
    print(f"Rotation: {f_rot:.4f} Hz  (period {rot_period:.3f} s)")

    print("Loading IQ...")
    x = load_iq(args.file, args.format, args.rate, args.max_seconds)
    duration = x.size / args.rate
    revs = duration * f_rot
    print(f"Loaded {x.size:,} samples = {duration:.1f} s = {revs:.1f} revolutions")
    if revs < 3:
        print("WARNING: fewer than ~3 revolutions recorded - the sine fit will be "
              "unreliable. Record longer (20+ revolutions ideal).")

    if args.peak_offset is None:
        peak_offset = find_peak_offset(x, args.rate)
        print(f"Strongest RFI peak at offset {peak_offset/1e3:.1f} kHz "
              f"({(args.center+peak_offset)/1e6:.4f} MHz)")
    else:
        peak_offset = args.peak_offset
        print(f"Tracking user peak at offset {peak_offset/1e3:.1f} kHz")

    print("Channelizing (shift peak to baseband + decimate)...")
    xc, crate = channelize(x, args.rate, peak_offset, args.channel_rate)
    print(f"Channel rate {crate:.0f} Hz, {xc.size:,} samples")

    slice_sec = rot_period / args.slices_per_rev
    print(f"Spectrogram slice = {slice_sec*1e3:.1f} ms "
          f"(freq resolution ~{1/slice_sec:.1f} Hz before interpolation)")
    f, t, Sxx, track = spectrogram_and_track(xc, crate, slice_sec, args.search_hz)

    fitres = fit_sine(t, track, f_rot, args.detrend_deg)
    fr, Y = track_periodogram(t, track, args.detrend_deg)

    # significance of the line at f_rot in the track periodogram
    k_rot = int(np.argmin(np.abs(fr - f_rot)))
    band = (fr > 0)
    prom = Y[k_rot] / (np.median(Y[band]) + 1e-20)

    # expected Doppler amplitude from geometry
    v_tip = 2 * np.pi * f_rot * args.radius
    expected_amp = v_tip / C * args.center

    print("\n================ RESULT ================")
    print(f"Fitted Doppler amplitude : {fitres['amplitude_hz']:.1f} Hz")
    print(f"Expected from geometry   : {expected_amp:.1f} Hz "
          f"(tip speed {v_tip:.1f} m/s at r={args.radius} m)")
    print(f"Sine-fit R^2             : {fitres['r2']:.3f}")
    print(f"Residual scatter         : {fitres['resid_rms_hz']:.1f} Hz")
    print(f"Periodogram prominence   : {prom:.1f}x median at the rotation freq")
    verdict = (fitres["r2"] > 0.5 and prom > 5)
    print("Verdict                  : " + (
        "DOPPLER SIGNATURE PRESENT - a sinusoid locks to the rotation frequency."
        if verdict else
        "No clear sinusoid at the rotation frequency (see plot; try a longer "
        "recording, a different peak, or check the rotation rate)."))
    print("========================================\n")

    make_plot(args.plot, f, t, Sxx, track, fitres, fr, Y, f_rot,
              args.center, peak_offset)


# ---------------------------------------------------------------------------
def _add_analysis_args(p):
    p.add_argument("file", help="raw IQ recording (taken while rotating)")
    p.add_argument("--rate", type=float, required=True, help="sample rate, Hz")
    p.add_argument("--center", type=float, required=True, help="center freq, Hz")
    p.add_argument("--rpm", type=float, help="ARM output rotation speed, RPM")
    p.add_argument("--rot-freq", type=float, help="rotation frequency, Hz (instead of --rpm)")
    p.add_argument("--radius", type=float, default=1.0, help="rotation radius, m")
    p.add_argument("--format", default="cf32", choices=["cf32", "u8", "i16"])
    p.add_argument("--max-seconds", type=float, default=30.0,
                   help="cap loaded duration (memory). Lower this on a small-RAM Pi.")
    p.add_argument("--peak-offset", type=float, default=None,
                   help="track this offset (Hz from center) instead of the strongest peak")
    p.add_argument("--channel-rate", type=float, default=4000.0,
                   help="decimated channel rate, Hz (covers +/- half this around the peak)")
    p.add_argument("--search-hz", type=float, default=800.0,
                   help="+/- window for tracking the peak around its center")
    p.add_argument("--slices-per-rev", type=float, default=20.0,
                   help="time slices per revolution (sets time resolution)")
    p.add_argument("--detrend-deg", type=int, default=2,
                   help="polynomial degree removed to suppress slow drift")
    p.add_argument("--plot", default="doppler_waterfall.png")


def main():
    ap = argparse.ArgumentParser(
        description="Rotational-Doppler: record raw IQ and/or analyze it for the "
                    "sine-at-rotation-frequency signature.")
    sub = ap.add_subparsers(dest="mode", required=True)

    # ---- record ----
    pr = sub.add_parser("record", help="capture raw IQ from the RTL-SDR to a file")
    pr.add_argument("file", help="output file (raw complex64 IQ)")
    pr.add_argument("--rate", type=float, default=2.4e6,
                    help="sample rate, Hz. Use 250000 + center on your peak for a small Pi.")
    pr.add_argument("--center", type=float, required=True, help="center freq, Hz")
    pr.add_argument("--seconds", type=float, required=True, help="record duration, s")
    pr.add_argument("--gain", default="40", help="tuner gain in dB, or 'auto'")
    pr.add_argument("--warmup", type=float, default=60.0,
                    help="discard this many seconds first so drift settles (default 60)")
    # optional: analyze right after recording
    pr.add_argument("--analyze", action="store_true",
                    help="analyze the recording immediately after capture")
    pr.add_argument("--rpm", type=float, help="(for --analyze) ARM rotation speed, RPM")
    pr.add_argument("--rot-freq", type=float, help="(for --analyze) rotation freq, Hz")
    pr.add_argument("--radius", type=float, default=1.0, help="(for --analyze) radius, m")

    # ---- analyze ----
    pa = sub.add_parser("analyze", help="analyze an existing raw IQ recording")
    _add_analysis_args(pa)

    args = ap.parse_args()

    if args.mode == "record":
        record_iq(args.file, args.rate, args.center, args.gain,
                  args.seconds, args.warmup)
        if args.analyze:
            # fill in the analysis defaults the record parser doesn't carry
            for name, default in (("format", "cf32"), ("max_seconds", args.seconds),
                                  ("peak_offset", None), ("channel_rate", 4000.0),
                                  ("search_hz", 800.0), ("slices_per_rev", 20.0),
                                  ("detrend_deg", 2), ("plot", "doppler_waterfall.png")):
                setattr(args, name, getattr(args, name, default))
            if not (args.rpm or args.rot_freq):
                sys.exit("--analyze needs --rpm or --rot-freq to fit against.")
            print()
            run_analysis(args)
        else:
            rate_hint = args.rpm if args.rpm else 100
            print("\nNext, analyze it with e.g.:")
            print(f"  python3 {sys.argv[0]} analyze {args.file} "
                  f"--rate {args.rate:g} --center {args.center:g} --rpm <your_RPM>")
    else:
        run_analysis(args)


if __name__ == "__main__":
    main()