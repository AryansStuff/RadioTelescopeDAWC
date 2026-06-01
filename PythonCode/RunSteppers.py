# I'm just putting this comment here as a manifestation that this code works, and that I don't need to debug for another day, I can't take much more of this and I let this comment serve as a prayer to the stepper motor gods.

import RPi.GPIO as GPIO
import time

PUL = 17
DIR = 27
DELAY = 0.003

GPIO.setwarnings(False)
GPIO.setmode(GPIO.BCM)
GPIO.setup(PUL, GPIO.OUT)
GPIO.setup(DIR, GPIO.OUT)
GPIO.output(DIR, GPIO.HIGH)
time.sleep(0.1)

print("Spinning")
try:
    while True:
        GPIO.output(PUL, GPIO.HIGH)
        time.sleep(DELAY)
        GPIO.output(PUL, GPIO.LOW)
        time.sleep(DELAY)
except KeyboardInterrupt:
    GPIO.cleanup()
    print("\nStopped.")