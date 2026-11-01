import argparse
import csv
import math
import queue
import signal
import sys
import threading
import time

import numpy as np
from scipy.signal import resample_poly

try:
    from rtlsdr import RtlSdr
except Exception:
    RtlSdr = None

try:
    import minimalmodbus
except Exception:
    minimalmodbus = None

from pilot_adf4351 import adf4351_registers, AdfPilot


REG_CONTROL = 0x8106
REG_SPEED = 0x8110
REG_ACTUAL = 0x8206
CTRL_STOP = 0x0000
CTRL_FWD = 0x0001
CTRL_REV = 0x0002


class SdrSource:
    def __init__(self, center_hz, fs_raw, gain, ppm=0, device_index=0):
        if RtlSdr is None:
            raise RuntimeError("pyrtlsdr not installed")
        self.sdr = RtlSdr(device_index=device_index)
        self.sdr.sample_rate = fs_raw
        self.sdr.center_freq = center_hz
        self.sdr.freq_correction = ppm if ppm else 1
        self.sdr.gain = gain
        self.fs_raw = fs_raw

    def read(self, n):
        return self.sdr.read_samples(n)

    def close(self):
        try:
            self.sdr.close()
        except Exception:
            pass


class MotorController:
    def __init__(self, port, slave=1, baud=9600):
        if minimalmodbus is None:
            raise RuntimeError("minimalmodbus not installed")
        self.inst = minimalmodbus.Instrument(port, slave)
        self.inst.serial.baudrate = baud
        self.inst.serial.timeout = 0.3
        self.inst.mode = minimalmodbus.MODE_RTU

    def set_rpm(self, rpm):
        self.inst.write_register(REG_SPEED, int(round(abs(rpm))), functioncode=6)

    def start(self, reverse=False):
        self.inst.write_register(REG_CONTROL, CTRL_REV if reverse else CTRL_FWD,
                                 functioncode=6)

    def stop(self):
        self.inst.write_register(REG_CONTROL, CTRL_STOP, functioncode=6)

    def read_rpm(self):
        raw = self.inst.read_register(REG_ACTUAL, functioncode=3)
        if raw >= 0x8000:
            raw -= 0x10000
        return raw

    def close(self):
        try:
            self.inst.serial.close()
        except Exception:
            pass


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
        self.src_idx = np.where((self.freqs >= src_band[0]) &
                                (self.freqs <= src_band[1]))[0]
        self.pil_idx = np.where((self.freqs >= pil_band[0]) &
                                (self.freqs <= pil_band[1]))[0]
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
            return self.freqs[k], mag[k]
        a, b, c = mag[k - 1], mag[k], mag[k + 1]
        denom = a - 2 * b + c
        if denom == 0:
            f_pk, pk = self.freqs[k], b
        else:
            d = 0.5 * (a - c) / denom
            f_pk, pk = self.freqs[k] + d * self.df, b - 0.25 * (a - c) * d
        noise = np.median(mag[idx]) + 1e-12
        return f_pk, 20.0 * np.log10(pk / noise)


class Acquisition:
    def __init__(self, args):
        self.args = args
        self.decim = int(round(args.fs_raw / args.fs_bb))
        self.fs_bb = args.fs_raw / self.decim
        self.shift = args.source_offset_hz
        self.pilot_bb_nominal = args.pilot_offset_hz - args.source_offset_hz
        self.q = queue.Queue(maxsize=64)
        self.stop_flag = threading.Event()
        self.t0 = None
        self.angle = 0.0
        self.last_motor_t = None
        nwin = int(round(args.win * self.fs_bb))
        self.tracker = DualTracker(
            self.fs_bb, nwin, nwin // 2,
            (-args.span / 2, args.span / 2),
            (self.pilot_bb_nominal - args.span / 2,
             self.pilot_bb_nominal + args.span / 2),
        )
        self.iq_file = open(args.out + ".c64", "wb")
        self.csv_file = open(args.out + ".csv", "w", newline="")
        self.csv = csv.writer(self.csv_file)
        self.csv.writerow(["t_s", "angle_rad", "rpm", "src_hz", "pilot_hz",
                           "corrected_hz", "src_snr", "pilot_snr"])
        self.motor = None
        self.adf = None

    def _producer(self, sdr):
        block = int(self.args.block_s * self.args.fs_raw)
        phase = 0.0
        w = 2.0 * np.pi * self.shift / self.args.fs_raw
        while not self.stop_flag.is_set():
            x = sdr.read(block)
            ts = time.time()
            n = len(x)
            k = np.arange(n)
            mix = np.exp(-1j * (phase + w * k))
            phase = (phase + w * n) % (2.0 * np.pi)
            bb = resample_poly(x * mix, 1, self.decim).astype(np.complex64)
            try:
                self.q.put((ts, bb), timeout=1.0)
            except queue.Full:
                pass

    def _read_motor(self):
        if self.motor is None:
            return 0.0
        try:
            return float(self.motor.read_rpm())
        except Exception:
            return float("nan")

    def _consumer(self):
        while not (self.stop_flag.is_set() and self.q.empty()):
            try:
                ts, bb = self.q.get(timeout=0.5)
            except queue.Empty:
                continue
            if self.t0 is None:
                self.t0 = ts
                self.last_motor_t = ts
            self.iq_file.write(bb.tobytes())
            rpm = self._read_motor()
            dt = ts - self.last_motor_t
            self.last_motor_t = ts
            if not np.isnan(rpm):
                self.angle += rpm * (2.0 * np.pi / 60.0) * dt
            pts = self.tracker.push(bb)
            last = None
            for i, (sf, ss, pf, ps) in enumerate(pts):
                corrected = sf - (pf - self.pilot_bb_nominal)
                t_pt = (ts - self.t0) + (i + 0.5) * (self.tracker.hop / self.fs_bb)
                self.csv.writerow([f"{t_pt:.4f}", f"{self.angle:.5f}",
                                   f"{rpm:.2f}", f"{sf:.5f}", f"{pf:.5f}",
                                   f"{corrected:.5f}", f"{ss:.2f}", f"{ps:.2f}"])
                last = (corrected, ss, ps)
            self.csv_file.flush()
            if last:
                sys.stdout.write(
                    f"\rt={ts - self.t0:7.1f}s rpm={rpm:6.1f} "
                    f"ang={np.degrees(self.angle) % 360:6.1f} "
                    f"corr={last[0]:+8.4f}Hz src_snr={last[1]:4.1f} "
                    f"pil_snr={last[2]:4.1f}")
                sys.stdout.flush()

    def run(self):
        if not self.args.no_adf:
            self.adf = AdfPilot(self.args.spi_bus, self.args.spi_dev,
                                self.args.adf_le_gpio, self.args.adf_lock_gpio)
            regs, actual = self.adf.program(self.args.pilot_hz)
            print("ADF4351 programmed: pilot =", f"{actual/1e6:.4f} MHz",
                  "regs =", [hex(r) for r in regs])
            if self.adf.wait_lock(self.args.lock_timeout):
                print("ADF4351 PLL locked.")
            else:
                print("WARNING: ADF4351 lock not detected; continuing.")

        sdr = SdrSource(self.args.center_hz, self.args.fs_raw,
                        self.args.gain, self.args.ppm)
        if not self.args.no_motor:
            self.motor = MotorController(self.args.motor_port,
                                        self.args.motor_slave)
        if self.motor and not self.args.no_move:
            self.motor.set_rpm(self.args.rpm)
            time.sleep(0.2)

        prod = threading.Thread(target=self._producer, args=(sdr,), daemon=True)
        cons = threading.Thread(target=self._consumer, daemon=True)
        prod.start()
        cons.start()
        try:
            self._run_schedule()
        except KeyboardInterrupt:
            pass
        finally:
            self.stop_flag.set()
            if self.motor and not self.args.no_move:
                self.motor.stop()
            prod.join(timeout=2.0)
            cons.join(timeout=5.0)
            if self.motor:
                self.motor.stop()
                self.motor.close()
            sdr.close()
            if self.adf:
                self.adf.close()
            self.iq_file.close()
            self.csv_file.close()
            print("\nstopped. wrote", self.args.out + ".c64", "and",
                  self.args.out + ".csv")

    def _hold(self, seconds):
        t_end = time.time() + seconds
        while time.time() < t_end and not self.stop_flag.is_set():
            time.sleep(0.1)

    def _run_schedule(self):
        # One steady rotation in a single direction: preroll (stationary) ->
        # spin at constant rpm for the sweep -> postroll. The analysis fits at a
        # single rotation rate, so the rotation must not change speed or reverse
        # mid-recording. Keep the sweep short enough that the cable service loop
        # is not overwound (a few turns), or add a slip ring for long runs.
        self._hold(self.args.preroll)
        if self.motor and not self.args.no_move:
            self.motor.start(reverse=False)
        self._hold(self.args.sweep)
        if self.motor and not self.args.no_move:
            self.motor.stop()
        self._hold(self.args.postroll)
        self.stop_flag.set()


def self_test():
    regs, actual = adf4351_registers(1420.405751e6)
    print("register calc for 1420.405751 MHz:")
    print("  actual pilot =", f"{actual/1e6:.5f} MHz")
    print("  regs =", [hex(r) for r in regs])

    fs_bb = 8000.0
    nwin = int(2.5 * fs_bb)
    pilot_nom = 3000.0
    tr = DualTracker(fs_bb, nwin, nwin // 2, (-4, 4),
                     (pilot_nom - 4, pilot_nom + 4))
    dur = 40.0
    t = np.arange(int(dur * fs_bb)) / fs_bb
    drift = 2.0 * np.sin(2 * np.pi * 0.02 * t + 0.5)
    doppler = 0.186 * np.sin(2 * np.pi / 120 * t + 0.7)
    src = np.exp(1j * 2 * np.pi * np.cumsum(drift + doppler) / fs_bb)
    pil = np.exp(1j * 2 * np.pi * np.cumsum(pilot_nom + drift) / fs_bb)
    rng = np.random.default_rng(0)
    x = (src + pil + 0.05 * (rng.standard_normal(len(t)) +
                             1j * rng.standard_normal(len(t)))).astype(np.complex64)
    sfs, corr = [], []
    for i in range(0, len(x), 8000):
        for sf, ss, pf, ps in tr.push(x[i:i + 8000]):
            sfs.append(sf)
            corr.append(sf - (pf - pilot_nom))
    sfs, corr = np.array(sfs), np.array(corr)
    print(f"\ndual-track over {len(sfs)} windows:")
    print(f"  raw source scatter (with drift) : {np.std(sfs)*1e3:7.1f} mHz")
    print(f"  drift-corrected scatter         : {np.std(corr)*1e3:7.1f} mHz")
    print(f"  drift removed by pilot          : {np.std(sfs)/max(np.std(corr),1e-9):6.1f}x")


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--center-hz", type=float, default=1420.405751e6)
    p.add_argument("--fs-raw", type=float, default=1.024e6)
    p.add_argument("--fs-bb", type=float, default=8000.0)
    p.add_argument("--gain", default="auto")
    p.add_argument("--ppm", type=int, default=0)
    p.add_argument("--source-offset-hz", type=float, default=0.0)
    p.add_argument("--pilot-offset-hz", type=float, default=3000.0)
    p.add_argument("--span", type=float, default=8.0)
    p.add_argument("--win", type=float, default=2.5)
    p.add_argument("--block-s", type=float, default=0.5)
    p.add_argument("--pilot-hz", type=float, default=1420.408751e6)
    p.add_argument("--spi-bus", type=int, default=0)
    p.add_argument("--spi-dev", type=int, default=0)
    p.add_argument("--adf-le-gpio", type=int, default=25)
    p.add_argument("--adf-lock-gpio", type=int, default=24)
    p.add_argument("--lock-timeout", type=float, default=5.0)
    p.add_argument("--no-adf", action="store_true")
    p.add_argument("--rpm", type=float, default=0.5)
    p.add_argument("--preroll", type=float, default=25.0)
    p.add_argument("--sweep", type=float, default=120.0)
    p.add_argument("--postroll", type=float, default=25.0)
    p.add_argument("--motor-port", default="/dev/ttyUSB0")
    p.add_argument("--motor-slave", type=int, default=1)
    p.add_argument("--no-motor", action="store_true")
    p.add_argument("--no-move", action="store_true")
    p.add_argument("--out", default="capture")
    p.add_argument("--self-test", action="store_true")
    return p


def main():
    args = build_parser().parse_args()
    if args.self_test:
        self_test()
        return
    if args.gain != "auto":
        args.gain = float(args.gain)
    acq = Acquisition(args)
    signal.signal(signal.SIGINT, lambda *_: acq.stop_flag.set())
    acq.run()


if __name__ == "__main__":
    main()