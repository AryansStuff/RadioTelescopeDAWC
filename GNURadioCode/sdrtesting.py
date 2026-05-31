import numpy as np
from rtlsdr import RtlSdr
import sounddevice as sd
from scipy import signal

# -------- SETTINGS --------
F_station = 96.7e6      # KISS FM
F_offset  = 250000      # offset-tune away from DC spike, mix back in software
Fs        = 1.2e6       # SDR sample rate
N         = 128000      # samples per block (multiple of 512, divisible by 25)
dec1      = 5           # 1.2 MHz -> 240 kHz
dec2      = 5           # 240 kHz -> 48 kHz audio
GAIN      = 'auto'      # or a number like 40
VOLUME    = 0.3
# --------------------------

Fs_demod   = Fs / dec1          # 240000  
audio_rate = Fs_demod / dec2    # 48000

# Setup SDR
sdr = RtlSdr()
sdr.sample_rate = Fs
sdr.center_freq = F_station - F_offset   # tune off to the side
sdr.gain = GAIN

# 75 us de-emphasis filter (US standard), one-pole IIR at the demod rate
tau   = 75e-6
alpha = np.exp(-1.0 / (Fs_demod * tau))
b_de, a_de = [1 - alpha], [1, -alpha]
zi_de = np.zeros(1)

# State carried between blocks for glitch-free output
prev_sample = 0.0 + 0.0j
lo_phase = 0.0
omega = -2 * np.pi * F_offset / Fs   # software local oscillator step

stream = sd.OutputStream(samplerate=int(audio_rate), channels=1, dtype='float32')
stream.start()

print(f"Tuned to {F_station/1e6:.1f} MHz — Ctrl+C to stop.")
try:
    while True:
        samples = sdr.read_samples(N)

        # Mix the station down to baseband (phase kept continuous across blocks)
        n = np.arange(len(samples))
        lo = np.exp(1j * (lo_phase + omega * n))
        lo_phase = (lo_phase + omega * len(samples)) % (2 * np.pi)
        x = samples * lo

        # Stage 1: low-pass + decimate to 240 kHz
        x = signal.decimate(x, dec1, ftype='fir')

        # FM demodulate (instantaneous frequency = angle of sample-to-sample product)
        x = np.concatenate(([prev_sample], x))
        prev_sample = x[-1]
        demod = np.angle(x[1:] * np.conj(x[:-1]))

        # De-emphasis
        demod, zi_de = signal.lfilter(b_de, a_de, demod, zi=zi_de)

        # Stage 2: decimate to 48 kHz audio
        audio = signal.decimate(demod, dec2, ftype='fir')

        # Scale, clip, play
        audio = np.clip(audio * VOLUME, -1.0, 1.0).astype(np.float32)
        stream.write(audio)
except KeyboardInterrupt:
    print("\nStopping.")
finally:
    stream.stop()
    stream.close()
    sdr.close()