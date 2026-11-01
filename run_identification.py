"""
Master run script for one DARFIM identification capture.

A single capture needs three things happening together:
  1. the ADF4351 pilot programmed and locked to the GPSDO,
  2. the central motor turning the arm at the target speed, and
  3. the horn motor counter-rotating so the horn stays on target,
while the SDR records the decimated baseband and the log.

The capture script (GNURadioCode/dopplerwaterfall.py) handles the pilot, the
central motor (over Modbus), the SDR recording, and the log. This script runs
the horn counter-rotation (PythonCode/horn_sync.py) alongside it on the same
preroll/sweep schedule, and stops both together.

    python3 run_identification.py --rpm 4 --out captures/run_4rpm
    python3 run_identification.py --rpm 0 --out captures/run_static   # the null

Flags after a bare -- are passed through to the capture script, e.g.
    python3 run_identification.py --rpm 4 --out run4 -- --gain 40 --ppm 1

Hardware only. Off the Pi the subprocesses exit immediately because their device
libraries are missing; that is expected.
"""

import argparse
import os
import signal
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CAPTURE = os.path.join(HERE, "GNURadioCode", "dopplerwaterfall.py")
HORN = os.path.join(HERE, "PythonCode", "horn_sync.py")


def main():
    ap = argparse.ArgumentParser(description="Run one identification capture")
    ap.add_argument("--rpm", type=float, required=True)
    ap.add_argument("--out", required=True, help="output tag (.c64 + .csv)")
    ap.add_argument("--preroll", type=float, default=25.0)
    ap.add_argument("--sweep", type=float, default=120.0)
    ap.add_argument("--postroll", type=float, default=25.0)
    ap.add_argument("--no-horn", action="store_true",
                    help="skip horn counter-rotation (central motor only)")
    ap.add_argument("passthrough", nargs="*",
                    help="extra args for the capture script, after --")
    args = ap.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    sched = ["--preroll", str(args.preroll), "--sweep", str(args.sweep)]

    horn_proc = None
    if args.rpm != 0 and not args.no_horn:
        horn_proc = subprocess.Popen(
            [sys.executable, HORN, "--rpm", str(args.rpm)] + sched)

    capture_cmd = ([sys.executable, CAPTURE, "--rpm", str(args.rpm),
                    "--out", args.out] + sched
                   + ["--postroll", str(args.postroll)] + args.passthrough)
    print("capture:", " ".join(capture_cmd))
    capture_proc = subprocess.Popen(capture_cmd)

    try:
        capture_proc.wait()
    except KeyboardInterrupt:
        capture_proc.send_signal(signal.SIGINT)
        capture_proc.wait()
    finally:
        if horn_proc is not None:
            horn_proc.send_signal(signal.SIGINT)
            try:
                horn_proc.wait(timeout=5.0)
            except subprocess.TimeoutExpired:
                horn_proc.kill()
    print("done:", args.out + ".c64", args.out + ".csv")


if __name__ == "__main__":
    main()
