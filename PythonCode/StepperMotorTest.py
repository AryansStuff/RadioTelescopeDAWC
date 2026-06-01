import RPi.GPIO as GPIO
import time

# --- Pins (BCM numbering, matching the 17 & 27 from `pinout`) ---
PUL = 17    # wired to PUL- on the driver
DIR = 27    # wired to DIR- on the driver

# --- Test settings ---
STEPS = 200        # at full-step, 200 steps = one full revolution
STEP_DELAY = 0.05  # seconds per edge; raise this to slow it down further

GPIO.setwarnings(False)
GPIO.setmode(GPIO.BCM)
GPIO.setup(PUL, GPIO.OUT)
GPIO.setup(DIR, GPIO.OUT)

GPIO.output(DIR, GPIO.HIGH)   # direction; flip to LOW to reverse
time.sleep(0.1)

try:
    print("Stepping... watch the shaft for a tick on each step.")
    for i in range(STEPS):
        GPIO.output(PUL, GPIO.HIGH)
        time.sleep(STEP_DELAY)
        GPIO.output(PUL, GPIO.LOW)
        time.sleep(STEP_DELAY)
        print(f"step {i + 1}/{STEPS}")
    print("Done.")
except KeyboardInterrupt:
    print("Stopped.")
finally:
    GPIO.cleanup()   # releases the pins so the hold/hum stops