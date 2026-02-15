from rtlsdr import RtlSdr
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import firwin, lfilter, iirfilter, sosfilt, freqz

# ---- Parameters ----
center_freq = 1.42041e9  # 21 cm Hydrogen line
sample_rate = 1.024e6  # RTL-SDR sample rate
output_rate = 32e3  # After decimation
fft_size = 32768
iir_alpha = 1e-3  # Single-pole IIR smoothing factor

# ---- Initialize SDR ----
sdr = RtlSdr()
sdr.sample_rate = sample_rate
sdr.center_freq = center_freq
sdr.gain = 'auto'

# ---- Design Xlating FIR Filter (low-pass) ----
decimation = int(sample_rate / output_rate)
cutoff = output_rate / 2
num_taps = 64
fir_taps = firwin(num_taps, cutoff / (sample_rate / 2))  # normalized
print(f"Decimation: {decimation}, FIR taps: {num_taps}")


# ---- Single-pole IIR (low-pass smoothing) ----
def single_pole_iir(x, alpha=iir_alpha, y_prev=0):
    y = np.zeros_like(x)
    for i in range(len(x)):
        y_prev = alpha * x[i] + (1 - alpha) * y_prev
        y[i] = y_prev
    return y


# ---- Real-time plot ----
plt.ion()
fig, ax = plt.subplots()
line, = ax.plot(np.linspace(-output_rate / 2, output_rate / 2, fft_size),
                np.zeros(fft_size))
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Power (dB)")
ax.set_title("Hydrogen Line Spectrum")

try:
    while True:
        # 1. Read samples
        samples = sdr.read_samples(fft_size * decimation)

        # 2. Apply FIR and decimate
        filtered = lfilter(fir_taps, 1.0, samples)[::decimation]

        # 3. FFT
        fft_vals = np.fft.fftshift(np.fft.fft(filtered, n=fft_size))

        # 4. Magnitude squared (power)
        power = np.abs(fft_vals) ** 2

        # 5. Single-pole IIR smoothing
        power_smooth = single_pole_iir(power)

        # 6. Log10 for dB
        power_db = 10 * np.log10(power_smooth + 1e-12)

        # 7. Update plot
        line.set_ydata(power_db)
        fig.canvas.draw()
        fig.canvas.flush_events()
except KeyboardInterrupt:
    print("Exiting...")
finally:
    sdr.close()
