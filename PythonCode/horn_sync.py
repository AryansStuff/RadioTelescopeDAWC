"""
Horn counter-rotation for a capture run.

While the central motor turns the arm at omega during the sweep, the horn motor
must turn at -omega so the horn keeps pointing at the target (d(beta)/dt =
-omega, derived in the paper; the horn's encoder reads beta, its angle relative
to the arm). This holds the horn still during the preroll, then ramps its target
angle at -omega for the sweep, matching the capture schedule.

    python3 horn_sync.py --rpm 4 --preroll 25 --sweep 120

Or import run_horn_sync(stop_event, rpm, preroll, sweep, ...) and start it in a
thread from the master script.

The motor pins, encoder pins, and counts-per-rev below must match the wiring.
PID gains almost certainly need tuning on the bench before a real run.
"""

import argparse
import threading
import time

from MotorTools import MotorController


def make_motor(fwd, back, enc_a, enc_b, cpr, kp, ki, kd, max_speed):
    return MotorController(
        motor_pin_forward=fwd,
        motor_pin_backward=back,
        encoder_a=enc_a,
        encoder_b=enc_b,
        counts_per_rev=cpr,
        kp=kp, ki=ki, kd=kd,
        max_speed=max_speed,
    )


def run_horn_sync(stop_event, rpm, preroll, sweep, motor=None, dt=0.02, **mk):
    """Hold during preroll, then counter-rotate at -rpm for the sweep.

    Ends when the sweep time is up or stop_event is set, whichever comes first.
    """
    if motor is None:
        motor = make_motor(**mk)
    deg_per_s = -rpm * 6.0                      # rpm -> deg/s, opposite sign
    hold_angle = motor.get_current_angle()
    motor.set_target_angle(hold_angle)

    # preroll: arm is stationary, so hold the horn in place
    t_end = time.monotonic() + preroll
    while time.monotonic() < t_end and not stop_event.is_set():
        motor.update()
        time.sleep(dt)

    # sweep: arm turns at +omega, so ramp the horn target at -omega
    t0 = time.monotonic()
    try:
        while not stop_event.is_set():
            elapsed = time.monotonic() - t0
            if elapsed > sweep:
                break
            motor.set_target_angle(hold_angle + deg_per_s * elapsed)
            motor.update()
            time.sleep(dt)
    finally:
        motor.stop()


def main():
    p = argparse.ArgumentParser(description="Counter-rotate the horn during a run")
    p.add_argument("--rpm", type=float, required=True, help="arm speed; horn runs at -rpm")
    p.add_argument("--preroll", type=float, default=25.0, help="stationary hold, s")
    p.add_argument("--sweep", type=float, default=120.0, help="rotation time, s")
    p.add_argument("--fwd", type=int, default=17, help="H-bridge forward pin (BCM)")
    p.add_argument("--back", type=int, default=18, help="H-bridge backward pin (BCM)")
    p.add_argument("--enc-a", type=int, default=22)
    p.add_argument("--enc-b", type=int, default=23)
    p.add_argument("--cpr", type=int, default=2048, help="encoder counts per rev")
    p.add_argument("--kp", type=float, default=0.02)
    p.add_argument("--ki", type=float, default=0.0)
    p.add_argument("--kd", type=float, default=0.005)
    p.add_argument("--max-speed", type=float, default=1.0)
    args = p.parse_args()

    stop = threading.Event()
    motor = make_motor(args.fwd, args.back, args.enc_a, args.enc_b, args.cpr,
                       args.kp, args.ki, args.kd, args.max_speed)
    print(f"horn: hold {args.preroll:.0f}s, then {-args.rpm:.2f} rpm for "
          f"{args.sweep:.0f}s (Ctrl+C to stop)")
    t = threading.Thread(target=run_horn_sync,
                         args=(stop, args.rpm, args.preroll, args.sweep, motor),
                         daemon=True)
    t.start()
    try:
        t.join()
    except KeyboardInterrupt:
        stop.set()
        t.join(timeout=2.0)
    finally:
        motor.stop()
        print("horn stopped.")


if __name__ == "__main__":
    main()
