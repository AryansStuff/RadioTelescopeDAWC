"""
Capture simulator for testing the analysis without the hardware.

Writes a decimated baseband .c64 and its companion .csv in the same format the
real capture produces, so analyzecapture.py and analyzeall.py run on it
unchanged. It models the things that matter for the pilot method:
  - a shared SDR clock drift (warm-up transient plus slow wander), several Hz,
    carried by both tones and removed by the pilot subtraction
  - a small independent phase-noise wander on each tone, which does not cancel
    and sets the residual floor
  - finite source and pilot SNR, a residual DC term, and one internal spur
  - rotation-speed jitter
The pilot carries the shared drift plus a tiny GPSDO residual instead of Doppler.

    python3 genv4cap.py --rpm 4 --out v4_4rpm
    python3 genv4cap.py --rpm 0 --out dawc_static        # the null
"""

import argparse
import csv
import numpy as np
from scipy.signal import butter, filtfilt

C = 299_792_458.0


def phase_noise(n, rms_hz, fs, rng):
    fn = rng.standard_normal(n)
    b, a = butter(2, 0.5 / (fs / 2))        # slow wander, ~0.5 Hz cutoff
    fn = filtfilt(b, a, fn)
    fn *= rms_hz / (np.std(fn) + 1e-12)
    return 2 * np.pi * np.cumsum(fn) / fs


def track(iq, fs, win_s, src_c, pil_c, half=8.0, pad=8):
    nwin = int(win_s * fs)
    hop = nwin // 2
    w = np.hanning(nwin)
    nfft = nwin * pad
    freqs = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs))
    df = freqs[1] - freqs[0]
    band = lambda c: np.where((freqs >= c - half) & (freqs <= c + half))[0]
    si, pi = band(src_c), band(pil_c)

    def pk(mag, idx):
        k = idx[int(np.argmax(mag[idx]))]
        if k <= 0 or k >= len(mag) - 1:
            return freqs[k], mag[k]
        a, b, c = mag[k - 1], mag[k], mag[k + 1]
        d = a - 2 * b + c
        s = 0.5 * (a - c) / d if d != 0 else 0
        return freqs[k] + s * df, b - 0.25 * (a - c) * s

    tw, sf, ss, pf, ps = [], [], [], [], []
    for s in range(0, len(iq) - nwin, hop):
        seg = iq[s:s + nwin] * w
        mag = np.abs(np.fft.fftshift(np.fft.fft(seg, n=nfft)))
        f1, p1 = pk(mag, si)
        f2, p2 = pk(mag, pi)
        tw.append((s + nwin / 2) / fs)
        sf.append(f1); ss.append(20 * np.log10(p1 / (np.median(mag[si]) + 1e-12)))
        pf.append(f2); ps.append(20 * np.log10(p2 / (np.median(mag[pi]) + 1e-12)))
    return map(np.array, (tw, sf, ss, pf, ps))


def build(args, rng):
    fs, n = args.fs, int(args.dur * args.fs)
    t = np.arange(n) / fs

    warm = args.drift * 0.5 * np.exp(-t / (args.dur * 0.4))
    wander = (0.6 * np.sin(2 * np.pi * 0.004 * t + rng.uniform(0, 6.28))
              + 0.4 * np.sin(2 * np.pi * 0.011 * t + rng.uniform(0, 6.28)))
    drift = warm + args.drift * 0.5 * wander

    if args.rpm > 0:
        om0 = 2 * np.pi * args.rpm / 60.0
        om = om0 * (1 + 0.02 * np.sin(2 * np.pi * 0.03 * t + rng.uniform(0, 6.28)))
        theta = np.cumsum(om) / fs + 0.7
        doppler = (om * args.radius / C * args.f0) * np.sin(theta)
    else:
        om = np.zeros(n)
        doppler = np.zeros(n)

    gpsdo = np.cumsum(rng.standard_normal(n)) * 1e-6
    gpsdo -= gpsdo.mean()
    gpsdo *= 0.014 / (np.std(gpsdo) + 1e-12)          # ~14 mHz residual

    pn_s = phase_noise(n, args.pn, fs, rng)
    pn_p = phase_noise(n, args.pn, fs, rng)

    src = np.exp(1j * (2 * np.pi * np.cumsum(args.src_off + drift + doppler) / fs + pn_s))
    pil = np.exp(1j * (2 * np.pi * np.cumsum(args.pilot_off + drift + gpsdo) / fs + pn_p))

    nwin = int(2.0 * fs)
    n_src = np.sqrt(nwin / 10 ** (args.snr_src / 10) / 2)
    a_src = 1.0
    a_pil = a_src * 10 ** ((args.snr_pilot - args.snr_src) / 20)
    noise = n_src * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    dc = args.dc * a_src * np.exp(1j * 0.4) * np.ones(n)
    spur = 0
    if args.spur is not None:
        spur = args.spur_lvl * a_src * np.exp(1j * (2 * np.pi * args.spur * t + 0.9))

    x = (a_src * src + a_pil * pil + dc + spur + noise).astype(np.complex64)
    return x, om


def main():
    p = argparse.ArgumentParser(description="Simulate a DARFIM capture for testing")
    p.add_argument("--rpm", type=float, default=4.0, help="0 = stationary null")
    p.add_argument("--fs", type=float, default=8000.0, help="baseband sample rate, Hz")
    p.add_argument("--f0", type=float, default=1420.406e6)
    p.add_argument("--radius", type=float, default=1.0)
    p.add_argument("--dur", type=float, default=30.0)
    p.add_argument("--src-off", type=float, default=0.0, help="signal baseband freq, Hz")
    p.add_argument("--pilot-off", type=float, default=3000.0, help="pilot baseband freq, Hz")
    p.add_argument("--snr-src", type=float, default=16.0)
    p.add_argument("--snr-pilot", type=float, default=19.0)
    p.add_argument("--drift", type=float, default=5.0, help="peak shared clock drift, Hz")
    p.add_argument("--pn", type=float, default=0.15, help="per-tone phase noise RMS, Hz")
    p.add_argument("--dc", type=float, default=0.12)
    p.add_argument("--spur", type=float, default=1800.0, help="spur location Hz; -1 disables")
    p.add_argument("--spur-lvl", type=float, default=0.10)
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--out", default=None, help="output tag; default v4_<rpm>rpm")
    args = p.parse_args()
    if args.spur is not None and args.spur < 0:
        args.spur = None

    rng = np.random.default_rng(args.seed)
    tag = args.out or f"v4_{int(args.rpm)}rpm"
    x, om = build(args, rng)
    x.tofile(f"{tag}.c64")

    n = len(x)
    win = 2.0 if args.rpm == 0 else float(np.clip((60.0 / args.rpm) / 12.0, 0.08, 4.0))
    a_th = (2 * np.pi * args.rpm / 60.0) * args.radius / C * args.f0 if args.rpm > 0 else 0.0
    tw, sf, ss, pf, ps = track(x, args.fs, win, args.src_off, args.pilot_off,
                               half=max(8.0, a_th + 4.0))
    corr = sf - args.src_off - (pf - args.pilot_off)
    idx = np.clip((tw * args.fs).astype(int), 0, n - 1)
    rpm_true = om[idx] * 60 / (2 * np.pi)
    rng2 = np.random.default_rng(args.seed + 7)
    rpm_log = np.round(rpm_true * 30) / 30 + rng2.normal(0, 0.004, len(rpm_true))
    dt = np.empty_like(tw)
    dt[0] = tw[0]
    dt[1:] = np.diff(tw)
    ang = np.cumsum(rpm_log * (2 * np.pi / 60) * dt)

    with open(f"{tag}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_s", "angle_rad", "rpm", "src_hz", "pilot_hz",
                    "corrected_hz", "src_snr", "pilot_snr"])
        for i in range(len(tw)):
            w.writerow([f"{tw[i]:.4f}", f"{ang[i]:.5f}", f"{rpm_log[i]:.3f}",
                        f"{sf[i]:.5f}", f"{pf[i]:.5f}", f"{corr[i]:.5f}",
                        f"{ss[i]:.2f}", f"{ps[i]:.2f}"])

    print(f"wrote {tag}.c64 ({x.nbytes/1e6:.2f} MB) and {tag}.csv "
          f"({n} samples, {args.dur:.0f}s, {args.rpm} rpm, fs {args.fs:.0f})")
    if args.rpm > 0:
        print(f"  theoretical Doppler amplitude: {a_th*1e3:.1f} mHz")


if __name__ == "__main__":
    main()
