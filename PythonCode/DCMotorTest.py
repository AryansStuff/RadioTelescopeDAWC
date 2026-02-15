import RPi.GPIO as GPIO
import time

IN1 = 22
IN2 = 23
ENA = 12  # PWM pin

GPIO.setmode(GPIO.BCM)
GPIO.setup(IN1, GPIO.OUT)
GPIO.setup(IN2, GPIO.OUT)
GPIO.setup(ENA, GPIO.OUT)

pwm = GPIO.PWM(ENA, 1000)  # 1 kHz PWM
pwm.start(0)

def forward(power):
    GPIO.output(IN1, 1)
    GPIO.output(IN2, 0)
    pwm.ChangeDutyCycle(power)

def reverse(power):
    GPIO.output(IN1, 0)
    GPIO.output(IN2, 1)
    pwm.ChangeDutyCycle(power)

def stop():
    pwm.ChangeDutyCycle(0)

try:
    print("Forward ramp")
    for power in range(0, 101, 10):
        print(f"Power: {power}%")
        forward(power)
        time.sleep(1)

    stop()
    time.sleep(2)

    print("Reverse ramp")
    for power in range(0, 101, 10):
        print(f"Power: {power}%")
        reverse(power)
        time.sleep(1)

    stop()

except KeyboardInterrupt:
    pass

pwm.stop()
GPIO.cleanup()
