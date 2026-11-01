"""
ADF4351 pilot-tone control over SPI from the Raspberry Pi.

The pilot is an ADF4351 synthesizer locked to the 10 MHz GPSDO reference. It is
injected through the directional coupler and used as the zero-Doppler frequency
reference that the analysis subtracts to remove the SDR clock drift.

This module does two things:
  1. adf4351_registers(fout, fref) computes the six 32-bit registers for a
     requested output frequency.
  2. AdfPilot writes them over SPI (LE on a GPIO) and reads the digital-lock
     pin so a run can confirm the PLL locked before recording.

Run it directly to program the pilot once and report lock:
    python3 pilot_adf4351.py --fout 1420.408751e6
    python3 pilot_adf4351.py --calc-only --fout 1420.408751e6   # no hardware

Register layout follows the ADF4351 datasheet. Confirm the computed value and
the lock state on the bench before trusting a capture.
"""

import argparse
import math
import time

try:
    import spidev
except Exception:
    spidev = None

try:
    import RPi.GPIO as GPIO
except Exception:
    GPIO = None


def adf4351_registers(fout, fref=10e6, chan_spacing=10e3, rf_power=3,
                      mux_lock=True):
    """Return (regs, actual_fout) for a requested output frequency.

    regs is R0..R5 as six integers. actual_fout is the frequency the PLL will
    actually produce, which can differ from fout by up to the channel spacing.
    """
    # pick the output divider so the VCO stays in its 2.2-4.4 GHz range
    rf_div_sel = 0
    div = 1
    while fout * div < 2.2e9 and div < 64:
        div *= 2
        rf_div_sel += 1
    vco = fout * div

    r_counter = 1
    fpfd = fref / r_counter
    mod = int(round(fpfd / chan_spacing))
    mod = max(2, min(mod, 4095))

    n = vco / fpfd
    int_n = int(n)
    frac = int(round((n - int_n) * mod)) % mod

    prescaler = 1 if int_n >= 75 else 0
    band_div = max(1, min(255, int(math.ceil(fpfd / 125e3))))
    muxout = 0b110 if mux_lock else 0b000     # 0b110 = digital lock detect
    charge_pump = 7

    r0 = (int_n << 15) | (frac << 3) | 0
    r1 = (prescaler << 27) | (1 << 15) | (mod << 3) | 1
    r2 = ((muxout << 26) | (r_counter << 14) | (charge_pump << 9) |
          (1 << 6) | 2)
    r3 = 0x000004B3
    r4 = ((1 << 23) | (rf_div_sel << 20) | (band_div << 12) | (1 << 5) |
          (rf_power << 3) | 4)
    r5 = 0x00580005

    actual = (int_n + frac / mod) * fpfd / div
    return [r0, r1, r2, r3, r4, r5], actual


class AdfPilot:
    """SPI driver for the ADF4351 board. LE and lock are plain GPIO lines."""

    def __init__(self, spi_bus=0, spi_dev=0, le_gpio=25, lock_gpio=24,
                 spi_hz=1_000_000):
        if spidev is None or GPIO is None:
            raise RuntimeError("spidev / RPi.GPIO not available (run on the Pi)")
        self.spi = spidev.SpiDev()
        self.spi.open(spi_bus, spi_dev)
        self.spi.max_speed_hz = spi_hz
        self.spi.mode = 0
        self.le = le_gpio
        self.lock = lock_gpio
        GPIO.setmode(GPIO.BCM)
        GPIO.setup(self.le, GPIO.OUT, initial=GPIO.HIGH)
        GPIO.setup(self.lock, GPIO.IN)

    def _write_reg(self, value):
        b = [(value >> 24) & 0xFF, (value >> 16) & 0xFF,
             (value >> 8) & 0xFF, value & 0xFF]
        GPIO.output(self.le, GPIO.LOW)
        self.spi.writebytes(b)
        GPIO.output(self.le, GPIO.HIGH)

    def program(self, fout, fref=10e6):
        """Write all six registers, high to low as the datasheet requires."""
        regs, actual = adf4351_registers(fout, fref)
        for value in reversed(regs):
            self._write_reg(value)
            time.sleep(0.001)
        return regs, actual

    def wait_lock(self, timeout=5.0):
        t_end = time.time() + timeout
        while time.time() < t_end:
            if GPIO.input(self.lock):
                return True
            time.sleep(0.05)
        return False

    def close(self):
        try:
            self.spi.close()
        except Exception:
            pass


def main():
    p = argparse.ArgumentParser(description="Program the ADF4351 pilot tone")
    p.add_argument("--fout", type=float, default=1420.408751e6,
                   help="pilot output frequency in Hz (default ~3 kHz above the line)")
    p.add_argument("--fref", type=float, default=10e6, help="GPSDO reference, Hz")
    p.add_argument("--spi-bus", type=int, default=0)
    p.add_argument("--spi-dev", type=int, default=0)
    p.add_argument("--le-gpio", type=int, default=25)
    p.add_argument("--lock-gpio", type=int, default=24)
    p.add_argument("--lock-timeout", type=float, default=5.0)
    p.add_argument("--calc-only", action="store_true",
                   help="print the register values and exit, no SPI")
    args = p.parse_args()

    regs, actual = adf4351_registers(args.fout, args.fref)
    print("requested  :", f"{args.fout/1e6:.6f} MHz")
    print("actual     :", f"{actual/1e6:.6f} MHz",
          f"(error {(actual-args.fout):+.1f} Hz)")
    print("registers  :", [hex(r) for r in regs])
    if args.calc_only:
        return

    adf = AdfPilot(args.spi_bus, args.spi_dev, args.le_gpio, args.lock_gpio)
    adf.program(args.fout, args.fref)
    print("locked     :", adf.wait_lock(args.lock_timeout))
    adf.close()


if __name__ == "__main__":
    main()
