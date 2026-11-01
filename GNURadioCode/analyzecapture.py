import argparse
import csv
import numpy as np


C = 299_792_458.0


class DualTracker:
    def __init__(self, fs_bb, nwin, hop, src_band, pil_band, pad=8):
        self.fs_bb = fs_bb
        self.nwin = nwin
        self.hop = hop
        self.win = np.hanning(nwin)
        self.win_gain = self.win.sum()
        self.nfft = int(nwin * pad)
        f = np.fft.fftfreq(self.nfft, d=1.0 / fs_bb)
        self.order = np.argsort(f)
        self.freqs = f[self.order]
        self.df = self.freqs[1] - self.freqs[0]
        self.src_idx = np.where((self.freqs >= src_band[0]) & (self.freqs <= src_band[1]))[0]
        self.pil_idx = np.where((self.freqs >= pil_band[0]) & (self.freqs <= pil_band[1]))[0]
        self.buf = np.zeros(0, dtype=np.complex64)

    def push(self, samples):
        self.buf = np.concatenate([self.buf, samples])
        out = []
        while len(self.buf) >= self.nwin:
            seg = self.buf[:self.nwin] * self.win
            mag = np.abs(np.fft.fft(seg, n=self.nfft) / self.win_gain)[self.order]
            sf, ss = self._peak(mag, self.src_idx)
            pf, ps = self._peak(mag, self.pil_idx)
            out.append((sf, ss, pf, ps))
            self.buf = self.buf[self.hop:]
        return out

    def _peak(self, mag, idx):
        k = idx[int(np.argmax(mag[idx]))]
        if k <= 0 or k >= len(mag) - 1:
            f_pk, pk = self.freqs[k], mag[k]
        else:
            a, b, c = mag[k - 1], mag[k], mag[k + 1]
            denom = a - 2 * b + c
            if denom == 0:
                f_pk, pk = self.freqs[k], b
            else:
                d = 0.5 * (a - c) / denom
                f_pk, pk = self.freqs[k] + d * self.df, b - 0.25 * (a - c) * d
        noise = np.median(mag[idx]) + 1e-12
        return f_pk, 20.0 * np.log10(pk / noise)


def load_iq(path):
    return np.fromfile(path, dtype=np.complex64)


def load_angle_csv(path):
    t, ang, rpm = [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            t.append(float(row["t_s"]))
            ang.append(float(row["angle_rad"]))
            rpm.append(float(row["rpm"]))
    return np.array(t), np.array(ang), np.array(rpm)


def track(iq, fs, win, span, pilot_off):
    nwin = int(round(win * fs))
    tr = DualTracker(fs, nwin, nwin // 2, (-span / 2, span / 2),
                     (pilot_off - span / 2, pilot_off + span / 2))
    rows = []
    for i in range(0, len(iq), 40000):
        rows += tr.push(iq[i:i + 40000])
    hop = nwin // 2
    tw = np.array([(j * hop + nwin / 2) / fs for j in range(len(rows))])
    sf = np.array([q[0] for q in rows])
    ss = np.array([q[1] for q in rows])
    pf = np.array([q[2] for q in rows])
    ps = np.array([q[3] for q in rows])
    return tw, sf, ss, pf, ps, tr


def fit(tw, corr, theta):
    X = np.column_stack([np.ones_like(tw), tw, tw ** 2, np.sin(theta), np.cos(theta)])
    beta, *_ = np.linalg.lstsq(X, corr, rcond=None)
    resid = corr - X @ beta
    amp = np.hypot(beta[3], beta[4])
    phase = np.arctan2(beta[4], beta[3])
    dof = len(corr) - 5
    cov = (np.sum(resid ** 2) / dof) * np.linalg.inv(X.T @ X)
    amp_err = np.sqrt((beta[3] / amp) ** 2 * cov[3, 3] + (beta[4] / amp) ** 2 * cov[4, 4]
                      + 2 * (beta[3] / amp) * (beta[4] / amp) * cov[3, 4])
    trend = beta[0] + beta[1] * tw + beta[2] * tw ** 2
    return amp, amp_err, phase, trend, resid


def main():
    ap = argparse.ArgumentParser(description="Analyze a DAWC capture (.c64 + .csv)")
    ap.add_argument("--iq", required=True, help="raw complex64 baseband IQ file")
    ap.add_argument("--csv", help="companion CSV with t_s, angle_rad, rpm (for fitting angle)")
    ap.add_argument("--fs", type=float, default=8000.0)
    ap.add_argument("--pilot-offset", type=float, default=3000.0)
    ap.add_argument("--span", type=float, default=8.0)
    ap.add_argument("--win", type=float, default=2.0)
    ap.add_argument("--f0", type=float, default=1420.405751e6)
    ap.add_argument("--out", default="analysis")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    iq = load_iq(args.iq)
    tw, sf, ss, pf, ps, tr = track(iq, args.fs, args.win, args.span, args.pilot_offset)
    corrected = sf - (pf - args.pilot_offset)

    if args.csv:
        ct, cang, crpm = load_angle_csv(args.csv)
        theta = np.interp(tw, ct, cang)
        omega = (cang[-1] - cang[0]) / (ct[-1] - ct[0])
        rpm_mean = crpm.mean()
    else:
        raise SystemExit("provide --csv for the rotation angle, or extend to accept --period")

    amp, amp_err, phase, trend, resid = fit(tw, corrected, theta)
    raw_amp, *_ = fit(tw, sf, theta)
    v = amp * C / args.f0
    v_err = amp_err * C / args.f0
    r = v / omega
    r_err = v_err / omega

    print("=" * 60)
    print("DAWC CAPTURE ANALYSIS")
    print("=" * 60)
    print(f"  file                : {args.iq}")
    print(f"  samples / duration  : {len(iq)} / {len(iq)/args.fs:.1f} s @ {args.fs:.0f} Hz")
    print(f"  windows tracked     : {len(tw)}")
    print(f"  mean rotation speed : {rpm_mean:.3f} rpm  (omega {omega:.4f} rad/s)")
    print(f"  median source SNR   : {np.median(ss):.1f} dB")
    print(f"  median pilot SNR    : {np.median(ps):.1f} dB")
    print(f"  raw Doppler amp     : {raw_amp*1e3:.1f} mHz  (drift-contaminated)")
    print(f"  corrected amplitude : {amp*1e3:.1f} +/- {amp_err*1e3:.1f} mHz")
    print(f"  tangential velocity : {v*100:.2f} +/- {v_err*100:.2f} cm/s")
    print(f"  recovered radius    : {r:.3f} +/- {r_err:.3f} m")
    print(f"  residual RMS        : {np.std(resid)*1e3:.1f} mHz")

    if args.no_plots:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(tw, sf, ".", ms=3, color="#d62728", alpha=0.5, label="raw source")
    ax.plot(tw, pf - args.pilot_offset, ".", ms=3, color="#7f7f7f", alpha=0.5, label="pilot drift")
    ax.plot(tw, corrected, ".", ms=4, color="#1f77b4", label="corrected")
    ax.set_xlabel("time (s)"); ax.set_ylabel("baseband frequency (Hz)")
    ax.set_title("Drift cancellation"); ax.legend(fontsize=8)
    fig.savefig(f"{args.out}_cancellation.png", dpi=150, bbox_inches="tight")

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(8, 5), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1]})
    a0.plot(tw, (corrected - trend) * 1e3, ".", ms=5, color="#1f77b4", label="corrected data")
    a0.plot(tw, amp * np.sin(theta) * 1e3, "r-", lw=1.2, label="sinusoid fit")
    a0.set_ylabel("Doppler (mHz)"); a0.legend()
    a0.set_title(f"Recovered Doppler (amp {amp*1e3:.1f} mHz, r {r:.3f} m)")
    a1.plot(tw, resid * 1e3, ".", ms=4, color="#555"); a1.axhline(0, color="k", lw=0.5)
    a1.set_ylabel("resid (mHz)"); a1.set_xlabel("time (s)")
    fig.savefig(f"{args.out}_doppler_fit.png", dpi=150, bbox_inches="tight")
    print(f"  figures             : {args.out}_cancellation.png, {args.out}_doppler_fit.png")


if __name__ == "__main__":
    main()