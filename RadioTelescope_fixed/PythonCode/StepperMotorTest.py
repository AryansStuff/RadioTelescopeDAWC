import RPi.GPIO as GPIO
import time

PUL = 17

GPIO.setwarnings(False)
GPIO.setmode(GPIO.BCM)
GPIO.setup(PUL, GPIO.OUT)

print("Slow PUL toggle")

try:
    state = False
    while True:
        state = not state
        GPIO.output(PUL, state)
        if state:
            print("PUL HIGH  ->  PUL- should read ~3.3 V")
        else:
            print("PUL LOW   ->  PUL- should read ~0 V")
        time.sleep(1.5)
except KeyboardInterrupt:
    GPIO.cleanup()
    print("\nStopped, pin released.")