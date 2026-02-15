import RPi.GPIO as GPIO
import time

# GPIO pin setup (your wiring)
IN1 = 22  # Red
IN2 = 23  # Blue
IN3 = 27  # Black
IN4 = 4   # Green

pins = [IN1, IN2, IN3, IN4]

# Setup GPIO
GPIO.setmode(GPIO.BCM)
for pin in pins:
    GPIO.setup(pin, GPIO.OUT)
    GPIO.output(pin, 0)  # turn off initially

# Full-step sequence (4 steps)
sequence = [
    (1, 0, 1, 0),
    (0, 1, 1, 0),
    (0, 1, 0, 1),
    (1, 0, 0, 1)
]

# Function to move motor
def step_motor(steps, delay=0.05, clockwise=True):
    seq_len = len(sequence)
    for i in range(steps):
        if clockwise:
            step = sequence[i % seq_len]
        else:
            step = sequence[(seq_len - (i % seq_len)) % seq_len]  # reverse
        for pin, val in zip(pins, step):
            GPIO.output(pin, int(val))  # ensure integer
        time.sleep(delay)

# Example usage
try:
    while True:
        # Rotate 1 revolution clockwise
        step_motor(200, delay=0.05, clockwise=True)
        time.sleep(1)
        # Rotate 1 revolution counterclockwise
        step_motor(200, delay=0.05, clockwise=False)
        time.sleep(1)
except KeyboardInterrupt:
    pass
finally:
    # Turn off all pins
    for pin in pins:
        GPIO.output(pin, 0)
    GPIO.cleanup()
