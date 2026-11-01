"""
Batch analysis of a set of captures (0/2/4/10 rpm, pilot running).

For each rotating run it writes a waterfall, a drift-cancellation plot, and a
Doppler sinusoid fit. Across all runs it writes a measured-vs-rpm chart, the
R^2 of measured against theory, the speed-to-speed linearity ratios, and a
summary CSV. The 0 rpm run is treated as the null.

Fitting is done at the rotation frequency known from the motor speed, not from
the logged angle, so it does not depend on the angle log being perfect.

    python3 analyzeall.py 0:dawc_static.c64 2:dawc_2rpm.c64 4:v4_4rpm.c64 10:dawc_10rpm.c64
    python3 analyzeall.py --fs 8000 --src-off 0 --pilot-off 3000 0:null.c64 4:run4.c64

Needs numpy + matplotlib.
"""

import argparse
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 299_792_458.0


def load_iq(path):
    x = np.fromfile(path, dtype=np.complex64)
    if x.size == 0:
        raise SystemExit(f"'{path}' is empty or not complex64.")
    return x


def theoretical_amp(rpm, radius, f0):
    return (2 * np.pi * rpm / 60.0) * radius * f0 / C


def track(iq, fs, win_s, src_c, pil_c, half_span, pad=8):
    """Dual-peak tracking per window with parabolic sub-bin interpolation."""
    nwin = max(16, int(win_s * fs))
    hop = max(1, nwin // 2)
    w = np.hanning(nwin)
    nfft = nwin * pad
    freqs = np.fft.fftshift(np.fft.fftfreq(nfft, 1.0 / fs))
    df = freqs[1] - freqs[0]

    def band(center):
        return np.where((freqs >= center - half_span) & (freqs <= center + half_span))[0]

    si, pi = band(src_c), band(pil_c)
    if len(si) == 0 or len(pi) == 0:
        raise SystemExit("Search band empty - check fs / src-off / pilot-off.")

    def peak(mag, idx):
        k = idx[int(np.argmax(mag[idx]))]
        if k <= 0 or k >= len(mag) - 1:
            return freqs[k], mag[k]
        a, b, c = mag[k - 1], mag[k], mag[k + 1]
        den = a - 2 * b + c
        d = 0.5 * (a - c) / den if den != 0 else 0.0
        return freqs[k] + d * df, b - 0.25 * (a - c) * d

    tw, sf, ss, pf, ps = [], [], [], [], []
    for s in range(0, len(iq) - nwin, hop):
        seg = iq[s:s + nwin] * w
        mag = np.abs(np.fft.fftshift(np.fft.fft(seg, n=nfft)))
        f1, p1 = peak(mag, si)
        f2, p2 = peak(mag, pi)
        nf1 = np.median(mag[si]) + 1e-12
        nf2 = np.median(mag[pi]) + 1e-12
        tw.append((s + nwin / 2) / fs)
        sf.append(f1); ss.append(20 * np.log10(max(p1, 1e-12) / nf1))
        pf.append(f2); ps.append(20 * np.log10(max(p2, 1e-12) / nf2))
    return (np.array(tw), np.array(sf), np.array(ss), np.array(pf), np.array(ps))


def fit_sinusoid(tw, corr, f_rot):
    """corr(t) = poly2 + A sin(2 pi f_rot t + phi). Return amp, err, trend, resid, offset."""
    w = 2 * np.pi * f_rot
    X = np.column_stack([np.ones_like(tw), tw, tw ** 2, np.sin(w * tw), np.cos(w * tw)])
    beta, *_ = np.linalg.lstsq(X, corr, rcond=None)
    resid = corr - X @ beta
    amp = np.hypot(beta[3], beta[4])
    dof = max(1, len(corr) - 5)
    cov = (np.sum(resid ** 2) / dof) * np.linalg.inv(X.T @ X)
    amp_err = 0.0
    if amp > 0:
        amp_err = np.sqrt((beta[3] / amp) ** 2 * cov[3, 3]
                          + (beta[4] / amp) ** 2 * cov[4, 4]
                          + 2 * (beta[3] / amp) * (beta[4] / amp) * cov[3, 4])
    trend = beta[0] + beta[1] * tw + beta[2] * tw ** 2
    return amp, amp_err, trend, resid, beta[0], beta


def win_for(rpm):
    if rpm <= 0:
        return 2.0
    return float(np.clip((60.0 / rpm) / 12.0, 0.08, 4.0))


def waterfall_fig(iq, fs, path, rpm, src_off, pilot_off):
    nwin = max(256, int(min(0.5, max(0.05, 4.0 / max(1, rpm))) * fs))
    hop = nwin // 2
    w = np.hanning(nwin)
    cols, times = [], []
    for s in range(0, len(iq) - nwin, hop):
        cols.append(np.fft.fftshift(np.abs(np.fft.fft(iq[s:s + nwin] * w))))
        times.append((s + nwin / 2) / fs)
    if not cols:
        return
    S = 20 * np.log10(np.array(cols).T + 1e-9)
    f = np.fft.fftshift(np.fft.fftfreq(nwin, 1.0 / fs))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.imshow(S, aspect="auto", origin="lower",
              extent=[times[0], times[-1], f[0], f[-1]],
              cmap="magma", vmin=S.max() - 55, vmax=S.max())
    ax.set_ylim(min(src_off, 0) - 300, max(abs(src_off), abs(pilot_off)) + 300)
    ax.axhline(src_off, color="cyan", ls=":", lw=0.6)
    ax.axhline(pilot_off, color="lime", ls=":", lw=0.6)
    ax.text(times[0], src_off + 40, "signal", color="cyan", fontsize=8)
    ax.text(times[0], pilot_off + 40, "pilot", color="lime", fontsize=8)
    ax.set_xlabel("time (s)"); ax.set_ylabel("baseband frequency (Hz)")
    ax.set_title(f"Recorded spectrum - {rpm} rpm (signal + pilot)")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def parse_runs(items):
    runs = []
    for it in items:
        if ":" not in it:
            raise SystemExit(f"run '{it}' must be rpm:path, e.g. 4:run4.c64")
        rpm, path = it.split(":", 1)
        runs.append((float(rpm), path))
    return runs


def main():
    ap = argparse.ArgumentParser(description="Batch DARFIM capture analysis")
    ap.add_argument("runs", nargs="+", help="rpm:path pairs, e.g. 0:null.c64 4:run4.c64")
    ap.add_argument("--fs", type=float, default=8000.0)
    ap.add_argument("--src-off", type=float, default=0.0)
    ap.add_argument("--pilot-off", type=float, default=3000.0)
    ap.add_argument("--f0", type=float, default=1420.406e6)
    ap.add_argument("--radius", type=float, default=1.0)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    runs = parse_runs(args.runs)
    os.makedirs(args.outdir, exist_ok=True)

    summary = []
    print("=" * 70)
    print("DARFIM ANALYSIS  (radius %.3g m, f0 %.6g Hz)" % (args.radius, args.f0))
    print("=" * 70)

    for rpm, path in runs:
        if not os.path.exists(path):
            print(f"  !! {path} not found - skipping {rpm} rpm")
            continue
        iq = load_iq(path)
        a_th = theoretical_amp(rpm, args.radius, args.f0)
        half = max(8.0, a_th + 4.0)
        win = win_for(rpm)
        tw, sf, ss, pf, ps = track(iq, args.fs, win, args.src_off, args.pilot_off, half)
        corrected = sf - args.src_off - (pf - args.pilot_off)
        tag = f"{int(rpm)}rpm"
        waterfall_fig(iq, args.fs, f"{args.outdir}/fig_waterfall_{tag}.png",
                      rpm, args.src_off, args.pilot_off)

        if rpm > 0:
            f_rot = rpm / 60.0
            amp, amp_err, trend, resid, offset, beta = fit_sinusoid(tw, corrected, f_rot)
            v = amp * C / args.f0
            r_meas = v / (2 * np.pi * f_rot)

            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.plot(tw, sf - args.src_off, ".", ms=3, color="#d62728", alpha=0.5,
                    label="raw signal (drift+Doppler)")
            ax.plot(tw, pf - args.pilot_off, ".", ms=3, color="#7f7f7f", alpha=0.5,
                    label="pilot drift")
            ax.plot(tw, corrected, ".", ms=4, color="#1f77b4", label="corrected")
            ax.set_xlabel("time (s)"); ax.set_ylabel("baseband frequency (Hz)")
            ax.set_title(f"Drift cancellation - {rpm} rpm"); ax.legend(fontsize=8)
            fig.savefig(f"{args.outdir}/fig_cancellation_{tag}.png", dpi=150,
                        bbox_inches="tight")
            plt.close(fig)

            fig, (a0, a1) = plt.subplots(2, 1, figsize=(8, 5), sharex=True,
                                         gridspec_kw={"height_ratios": [3, 1]})
            a0.plot(tw, (corrected - trend) * 1e3, ".", ms=5, color="#1f77b4",
                    label="corrected data")
            tt = np.linspace(tw[0], tw[-1], 1000)
            w = 2 * np.pi * f_rot
            a0.plot(tt, (beta[3] * np.sin(w * tt) + beta[4] * np.cos(w * tt)) * 1e3,
                    "r-", lw=1.2, label="sinusoid fit")
            a0.set_ylabel("Doppler (mHz)"); a0.legend()
            a0.set_title(f"Recovered Doppler - {rpm} rpm  (amp {amp*1e3:.0f} mHz, r {r_meas:.3f} m)")
            a1.plot(tw, resid * 1e3, ".", ms=4, color="#555")
            a1.axhline(0, color="k", lw=0.5)
            a1.set_ylabel("resid (mHz)"); a1.set_xlabel("time (s)")
            fig.savefig(f"{args.outdir}/fig_fit_{tag}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

            print(f"\n  {rpm:>4g} rpm | win {win*1000:4.0f} ms | theory {a_th*1e3:7.1f} mHz | "
                  f"measured {amp*1e3:7.1f} +/- {amp_err*1e3:4.1f} mHz | "
                  f"err {100*(amp-a_th)/a_th:+5.1f}% | r {r_meas:.3f} m | "
                  f"offset {offset*1e3:+.1f} mHz")
            summary.append(dict(rpm=rpm, theory=a_th, measured=amp, err=amp_err,
                                r=r_meas, offset=offset, snr=np.median(ss)))
        else:
            p = np.polyfit(tw, corrected, 2)
            null_level = np.std(corrected - np.polyval(p, tw))
            print(f"\n  {rpm:>4g} rpm | NULL: residual scatter {null_level*1e3:.1f} mHz "
                  f"(no rotation-locked signal)")
            summary.append(dict(rpm=0, theory=0.0, measured=null_level, err=0.0,
                                r=float("nan"), offset=np.mean(corrected),
                                snr=np.median(ss)))

    rot = [s for s in summary if s["rpm"] > 0]
    if rot:
        rpm_a = np.array([s["rpm"] for s in rot])
        meas_a = np.array([s["measured"] for s in rot])
        th_a = np.array([s["theory"] for s in rot])
        ss_res = np.sum((meas_a - th_a) ** 2)
        ss_tot = np.sum((meas_a - np.mean(meas_a)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

        fig, ax = plt.subplots(figsize=(8, 5))
        xx = np.linspace(0, max(rpm_a) * 1.1, 50)
        ax.plot(xx, [theoretical_amp(x, args.radius, args.f0) * 1e3 for x in xx],
                "r--", lw=1.2, label="theory")
        ax.errorbar(rpm_a, meas_a * 1e3, yerr=[s["err"] * 1e3 for s in rot],
                    fmt="o", ms=8, capsize=4, color="#1f77b4", label="measured")
        null = [s for s in summary if s["rpm"] == 0]
        if null:
            ax.plot(0, null[0]["measured"] * 1e3, "s", ms=9, color="#2ca02c",
                    label="null (0 rpm)")
        ax.set_xlabel("rotational speed (rpm)"); ax.set_ylabel("Doppler shift (mHz)")
        ax.set_title("Measured vs. theoretical Doppler shift (R2 = %.4f)" % r2)
        ax.legend(); ax.grid(alpha=0.3)
        fig.savefig(f"{args.outdir}/fig_shift_vs_rpm.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        print("\n" + "-" * 70)
        print("  MEASURED vs THEORY:  R2 = %.4f" % r2)
        order = np.argsort(rpm_a)
        for i in range(1, len(order)):
            lo, hi = order[i - 1], order[i]
            print("  linearity  A(%g)/A(%g) = %.3f  (expect %.3f)" %
                  (rpm_a[hi], rpm_a[lo], meas_a[hi] / meas_a[lo], rpm_a[hi] / rpm_a[lo]))
        offs = np.array([s["offset"] for s in rot])
        print("  residual offset after cancellation: mean %.1f mHz (|max| %.1f mHz)"
              % (np.mean(offs) * 1e3, np.max(np.abs(offs)) * 1e3))

    with open(f"{args.outdir}/summary.csv", "w", newline="") as f:
        wtr = csv.writer(f)
        wtr.writerow(["rpm", "theoretical_hz", "measured_hz", "measured_err_hz",
                      "recovered_radius_m", "residual_offset_hz", "median_snr_db"])
        for s in summary:
            wtr.writerow([s["rpm"], f"{s['theory']:.4f}", f"{s['measured']:.4f}",
                          f"{s['err']:.4f}", f"{s['r']:.4f}", f"{s['offset']:.4f}",
                          f"{s['snr']:.2f}"])
    print("\n  figures + summary.csv written to", os.path.abspath(args.outdir))
    print("=" * 70)


if __name__ == "__main__":
    main()
