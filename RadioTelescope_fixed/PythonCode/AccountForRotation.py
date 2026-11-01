import time
from MotorTools import MotorController

# --- User input ---
rpm = float(input("What speed in RPM will this motor run? "))
time_int = float(input("What time interval should corrections happen between (s)? "))

# Convert RPM → degrees
delta_angle = rpm * 6.0 * time_int

# --- Motor setup ---
motor = MotorController(
    motor_pin_forward=17,
    motor_pin_backward=18,
    encoder_a=22,
    encoder_b=23,
    counts_per_rev=200,
    kp=0.02,
    kd=0.005
)

# Set relative motion
current_angle = motor.get_current_angle()
motor.set_target_angle(current_angle + delta_angle)

# --- Control loop ---
try:
    while True:
        error, control = motor.update()

        if abs(error) < 0.5:  # 0.5 degree tolerance
            break

        time.sleep(time_int)

finally:
    motor.stop()
