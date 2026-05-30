#!/usr/bin/env python3
"""
record_spectrum.py
Reproduces the GNU Radio data-analysis sequence from Figure 4 of the DAWC paper:

    Soapy RTL-SDR Source
      -> Frequency Xlating FIR Filter   (cut broadband noise / select the band)
      -> Stream to Vector (8192)
      -> FFT (8192)                      (frequency domain)
      -> Complex to Mag^2                (power)
      -> Single Pole IIR Filter          (time averaging, decaying memory)
      -> Log10 (x10)                     (scale to dB)
      -> File Sink                       (save spectra to disk)

Run it ONCE pointed at the sun while STATIONARY, and ONCE while ROTATING,
writing to two different files. Then analyse with analyze_shift.py.

    python3 record_spectrum.py stationary.dat --seconds 600
    python3 record_spectrum.py rotating.dat   --seconds 600

NOTE: GNU Radio block constructors vary slightly between versions (this targets
3.10). If the soapy source line errors, see the osmosdr alternative noted below.
"""

import time
import argparse

from gnuradio import gr, blocks, fft, filter
from gnuradio.fft import window
from gnuradio import soapy


class spectrum_recorder(gr.top_block):
    def __init__(self, outfile, samp_rate, center_freq, gain,
                 fft_size, alpha, bw, keep_one_in_n):
        gr.top_block.__init__(self, "DAWC Spectrum Recorder")

        # ---- Soapy RTL-SDR Source ("Take input from the RTL-SDR") ------
        dev = 'driver=rtlsdr'
        self.src = soapy.source(dev, "fc32", 1, '', '', [''], [''])
        self.src.set_sample_rate(0, samp_rate)
        self.src.set_gain_mode(0, False)            # manual gain (AGC off)
        self.src.set_frequency(0, center_freq)
        self.src.set_frequency_correction(0, 0)
        self.src.set_gain(0, 'TUNER', gain)
        # osmosdr alternative (if you don't have gr-soapy):
        #   import osmosdr
        #   self.src = osmosdr.source(args="rtl=0")
        #   self.src.set_sample_rate(samp_rate); self.src.set_center_freq(center_freq)
        #   self.src.set_gain(gain)

        # ---- Frequency Xlating FIR Filter ------------------------------
        # Low-pass taps select the band of interest and cut broadband noise
        # ("Cut out of broadband noise and remove uncertainty").
        taps = filter.firdes.low_pass(1.0, samp_rate, bw / 2.0,
                                      bw / 10.0, window.WIN_HAMMING)
        self.xlate = filter.freq_xlating_fir_filter_ccc(
            1, taps, 0.0, samp_rate)                # decim=1, no freq offset

        # ---- Stream to Vector (group samples into 8192-point frames) ---
        self.s2v = blocks.stream_to_vector(gr.sizeof_gr_complex, fft_size)

        # ---- FFT (to the frequency domain) -----------------------------
        self.fft = fft.fft_vcc(fft_size, True,
                               window.blackmanharris(fft_size), True, 1)

        # ---- Complex to Mag^2 ("...in terms of power") -----------------
        self.mag2 = blocks.complex_to_mag_squared(fft_size)

        # ---- Single Pole IIR Filter (FFT averaging over time) ----------
        # "Averages values from over time with past values having
        #  decreasing relevancy." Smaller alpha = smoother / slower.
        self.iir = filter.single_pole_iir_filter_ff(alpha, fft_size)

        # ---- Log10 (x10 -> dB, "scale it down") ------------------------
        self.log = blocks.nlog10_ff(10, fft_size, 0)

        # ---- Thin the frame rate so the file stays a sane size ---------
        # The IIR filter already does the averaging, so we only keep a few
        # converged spectra per second instead of every single frame.
        self.keep = blocks.keep_one_in_n(gr.sizeof_float * fft_size,
                                         keep_one_in_n)

        # ---- File Sink ("Save it to a file") ---------------------------
        self.sink = blocks.file_sink(gr.sizeof_float * fft_size, outfile, False)
        self.sink.set_unbuffered(False)

        # ---- Wire the chain together -----------------------------------
        self.connect(self.src, self.xlate, self.s2v, self.fft,
                     self.mag2, self.iir, self.log, self.keep, self.sink)


def main():
    p = argparse.ArgumentParser(description="Record averaged spectra (paper Figure 4).")
    p.add_argument('outfile', help='output file, e.g. stationary.dat or rotating.dat')
    p.add_argument('--seconds', type=float, default=600.0,
                   help='record time in seconds (default 600 = 10 min, as in the paper)')
    p.add_argument('--samp-rate', type=float, default=2.4e6)
    p.add_argument('--center-freq', type=float, default=1.42041e9,
                   help='hydrogen-line band centre (Hz)')
    p.add_argument('--gain', type=float, default=40.0)
    p.add_argument('--fft-size', type=int, default=8192,
                   help='8192 bins, as in the paper')
    p.add_argument('--alpha', type=float, default=0.01,
                   help='IIR averaging coefficient (smaller = smoother)')
    p.add_argument('--bw', type=float, default=2.0e6,
                   help='bandwidth kept by the xlating filter (Hz)')
    args = p.parse_args()

    frame_rate = args.samp_rate / args.fft_size
    keep_one_in_n = max(1, int(frame_rate / 5))     # ~5 spectra/sec to disk
    bin_width = args.samp_rate / args.fft_size

    tb = spectrum_recorder(args.outfile, args.samp_rate, args.center_freq,
                           args.gain, args.fft_size, args.alpha, args.bw,
                           keep_one_in_n)
    print(f"Recording {args.seconds:.0f}s -> {args.outfile}  "
          f"(bin width = {bin_width:.0f} Hz). Ctrl+C to stop early.")
    tb.start()
    try:
        time.sleep(args.seconds)
    except KeyboardInterrupt:
        pass
    tb.stop()
    tb.wait()
    print("Done.")


if __name__ == '__main__':
    main()