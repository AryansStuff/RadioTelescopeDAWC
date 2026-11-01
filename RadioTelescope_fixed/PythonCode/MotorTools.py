from gpiozero import Motor, RotaryEncoder
import time


class PID:
    def __init__(self, kp, ki, kd, output_limits=(-1.0, 1.0)):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.min_out, self.max_out = output_limits

        self.integral = 0.0
        self.prev_error = None

    def reset(self):
        self.integral = 0.0
        self.prev_error = None

    def compute(self, error, dt):
        if dt <= 0:
            return 0.0

        self.integral += error * dt

        derivative = 0.0
        if self.prev_error is not None:
            derivative = (error - self.prev_error) / dt

        self.prev_error = error

        output = (
            self.kp * error +
            self.ki * self.integral +
            self.kd * derivative
        )

        # Clamp output
        return max(self.min_out, min(self.max_out, output))


class MotorController:
    def __init__(
        self,
        motor_pin_forward,
        motor_pin_backward,
        encoder_a,
        encoder_b,
        counts_per_rev,
        kp=0.01,
        ki=0.0,
        kd=0.0,
        max_speed=1.0
    ):
        # Motor and encoder
        self.motor = Motor(motor_pin_forward, motor_pin_backward)
        self.encoder = RotaryEncoder(encoder_a, encoder_b, max_steps=0)

        # Encoder → angle
        self.counts_per_rev = counts_per_rev
        self.deg_per_count = 360.0 / counts_per_rev

        # PID controller
        self.pid = PID(kp, ki, kd, output_limits=(-max_speed, max_speed))

        # State
        self.target_angle = 0.0
        self.last_time = time.monotonic()

    # ----------------------------
    # Angle interface
    # ----------------------------
    def set_target_angle(self, angle_deg):
        """Set desired angle in degrees."""
        self.target_angle = angle_deg
        self.pid.reset()

    def get_current_angle(self):
        """Return current angle in degrees."""
        return self.encoder.steps * self.deg_per_count

    # ----------------------------
    # Control update
    # ----------------------------
    def update(self):
        """
        Call this repeatedly in your main loop.
        Returns (error, control_output)
        """
        now = time.monotonic()
        dt = now - self.last_time
        self.last_time = now

        current_angle = self.get_current_angle()
        error = self.target_angle - current_angle

        control = self.pid.compute(error, dt)
        self._apply_motor(control)

        return error, control

    # ----------------------------
    # Motor output
    # ----------------------------
    def _apply_motor(self, value):
        if value > 0:
            self.motor.forward(value)
        elif value < 0:
            self.motor.backward(-value)
        else:
            self.motor.stop()

    # ----------------------------
    # Safety / utility
    # ----------------------------
    def stop(self):
        self.motor.stop()
        self.pid.reset()


# -------------------------------------------------
# Example usage (can be deleted in production)
# -------------------------------------------------
if __name__ == "__main__":
    motor = MotorController(
        motor_pin_forward=17,
        motor_pin_backward=18,
        encoder_a=22,
        encoder_b=23,
        counts_per_rev=2048,
        kp=0.02,
        kd=0.005
    )

    motor.set_target_angle(45.0)

    try:
        while True:
            error, control = motor.update()
            print(f"Angle: {motor.get_current_angle():.2f}° | Error: {error:.2f} | Output: {control:.2f}")
            time.sleep(0.02)

    except KeyboardInterrupt:
        motor.stop()
        print("Stopped safely.")
