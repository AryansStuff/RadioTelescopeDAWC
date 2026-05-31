import RPi.GPIO as GPIO
import time

PUL = 17   # wired to PUL-
DIR = 27   # wired to DIR-

GPIO.setmode(GPIO.BCM)
GPIO.setup(PUL, GPIO.OUT)
GPIO.setup(DIR, GPIO.OUT)

GPIO.output(DIR, GPIO.HIGH)   # direction; flip to LOW to reverse

while True:                   # spin continuously
    GPIO.output(PUL, GPIO.HIGH)
    time.sleep(0.0005)        # smaller delay = faster
    GPIO.output(PUL, GPIO.LOW)
    time.sleep(0.0005)