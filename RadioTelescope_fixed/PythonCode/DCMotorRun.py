import RPi.GPIO as GPIO
import time

# L298N pins (BCM) — change to match your wiring
IN1 = 22   # direction A
IN2 = 23   # direction B
ENA = 12   # enable / speed pin

GPIO.setwarnings(False)
GPIO.setmode(GPIO.BCM)
for p in (IN1, IN2, ENA):
    GPIO.setup(p, GPIO.OUT)

# Set direction (one way)
GPIO.output(IN1, GPIO.HIGH)
GPIO.output(IN2, GPIO.LOW)
# Full power: ENA HIGH = 100% (no PWM needed)
GPIO.output(ENA, GPIO.HIGH)

print("Motor at full power. Ctrl+C to stop.")
try:
    while True:
        time.sleep(1)
except KeyboardInterrupt:
    GPIO.output(ENA, GPIO.LOW)
    GPIO.cleanup()
    print("\nStopped.")