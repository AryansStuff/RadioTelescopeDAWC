from rtlsdr import RtlSdr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# -------- SETTINGS --------
CENTER_FREQ = 1.42041e9
SAMPLE_RATE = 1.024e6
FFT_SIZE = 4096
GAIN = 30
# --------------------------

# Setup SDR
sdr = RtlSdr()
sdr.sample_rate = SAMPLE_RATE
sdr.center_freq = CENTER_FREQ
sdr.gain = GAIN

# Frequency axis
freqs = np.fft.fftshift(np.fft.fftfreq(FFT_SIZE, 1/SAMPLE_RATE))
absolute_freqs = freqs + CENTER_FREQ

# Setup plot
fig, ax = plt.subplots()
line, = ax.plot(absolute_freqs, np.zeros(FFT_SIZE))
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Power (dB)")
ax.set_title("Real-Time Spectrum")
ax.set_ylim(-80, 20)

def update(frame):
    samples = sdr.read_samples(FFT_SIZE)

    fft_vals = np.fft.fftshift(np.fft.fft(samples))
    power = np.abs(fft_vals)**2
    power_db = 10*np.log10(power + 1e-12)

    line.set_ydata(power_db)
    return line,

ani = FuncAnimation(fig, update, interval=50, blit=True)

try:
    plt.show()
except KeyboardInterrupt:
    pass
finally:
    sdr.close()
