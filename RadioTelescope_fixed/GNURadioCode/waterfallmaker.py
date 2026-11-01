import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_iq(path):
    x = np.fromfile(path, dtype=np.complex64)
    if x.size == 0:
        raise SystemExit(f"'{path}' is empty or not readable as complex64 IQ.")
    return x


def spectrogram(x, fs, win_s, hop_frac=0.5):
    nwin = int(win_s * fs)
    hop = max(1, int(nwin * hop_frac))
    w = np.hanning(nwin)
    cols, times = [], []
    for s in range(0, len(x) - nwin, hop):
        seg = x[s:s + nwin] * w
        cols.append(np.fft.fftshift(np.abs(np.fft.fft(seg))))
        times.append((s + nwin / 2) / fs)
    S = 20 * np.log10(np.array(cols).T + 1e-9)
    f = np.fft.fftshift(np.fft.fftfreq(nwin, 1 / fs))
    return f, np.array(times), S


def main():
    ap = argparse.ArgumentParser(description="Waterfall (spectrogram) from a .c64 capture")
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, default=8000.0)
    ap.add_argument("--win", type=float, default=0.5, help="display FFT window (s)")
    ap.add_argument("--pilot-offset", type=float, default=3000.0)
    ap.add_argument("--zoom", type=float, default=12.0, help="+/- Hz around source for zoom panel")
    ap.add_argument("--out", default="waterfall")
    a = ap.parse_args()

    x = load_iq(a.iq)
    f, t, S = spectrogram(x, a.fs, a.win)

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))
    a0.imshow(S, aspect="auto", origin="lower", extent=[t[0], t[-1], f[0], f[-1]],
              cmap="magma", vmin=S.max() - 55, vmax=S.max())
    a0.set_ylim(-a.zoom - 5, a.pilot_offset + 60)
    a0.axhline(0, color="cyan", ls=":", lw=0.6)
    a0.axhline(a.pilot_offset, color="lime", ls=":", lw=0.6)
    a0.text(t[0] + 0.5, 60, "source (carries Doppler)", color="cyan", fontsize=8)
    a0.text(t[0] + 0.5, a.pilot_offset + 90, "pilot (reference)", color="lime", fontsize=8)
    a0.set_title("Full baseband: two tones")
    a0.set_xlabel("time (s)"); a0.set_ylabel("frequency (Hz)")

    im = a1.imshow(S, aspect="auto", origin="lower", extent=[t[0], t[-1], f[0], f[-1]],
                   cmap="magma", vmin=S.max() - 40, vmax=S.max())
    a1.set_ylim(-a.zoom, a.zoom)
    a1.set_title("Zoom on the source tone: frequency vs time")
    a1.set_xlabel("time (s)"); a1.set_ylabel("frequency (Hz)")
    fig.colorbar(im, ax=a1, label="power (dB)")
    fig.tight_layout()
    fig.savefig(a.out + ".png", dpi=150, bbox_inches="tight")
    print(f"wrote {a.out}.png  ({len(t)} time slices, {a.win*1000:.0f} ms window, "
          f"{a.fs/int(a.win*a.fs):.1f} Hz resolution)")


if __name__ == "__main__":
    main()