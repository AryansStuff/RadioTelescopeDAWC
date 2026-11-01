# Radio Telescope - DARFIM identification pipeline

Spin the antenna slowly; a fixed-direction interferer gets a rotation-locked
Doppler wobble while the on-axis target does not. A GPSDO-locked pilot tone,
injected by cable, cancels the receiver's clock drift so the sub-hertz wobble
is recoverable at 2-10 rpm. This repo captures a run and recovers that
fingerprint. (Removing the RFI is a later step and is not in here.)

See `MANIFEST.md` for what every file is.

## 1. Setup (once)

```
# system packages + python deps
bash setup.sh                 # or: sudo apt install $(cat apt-packages.txt); pip install -r requirements.txt
```

On the Raspberry Pi also enable SPI (`sudo raspi-config` -> Interface -> SPI)
for the ADF4351 pilot.

Nothing in the capture/motor path runs off the Pi - the device libraries
(`pyrtlsdr`, `minimalmodbus`, `spidev`, `RPi.GPIO`, `gpiozero`) are only present
there. The analysis and the simulator run anywhere with numpy + scipy +
matplotlib.

## 2. Bench checks (optional but recommended)

```
python3 GNURadioCode/pilot_adf4351.py --calc-only --fout 1420.408751e6   # register math, no hardware
python3 GNURadioCode/pilot_adf4351.py --fout 1420.408751e6               # program + confirm lock (on Pi)
python3 GNURadioCode/rtlsdrtest.py                                       # see the band live
python3 PythonCode/StepperMotorTest.py                                   # central-motor wiring
python3 PythonCode/DCMotorTest.py                                        # horn-motor wiring
```

## 3. Take the data (on the Pi)

Point the rotation axis at the target (the Sun) first. Then take one run per
speed - the null plus a few low speeds:

```
python3 run_identification.py --rpm 0  --out captures/run_static
python3 run_identification.py --rpm 2  --out captures/run_2rpm
python3 run_identification.py --rpm 4  --out captures/run_4rpm
python3 run_identification.py --rpm 10 --out captures/run_10rpm
```

Each writes `captures/run_*.c64` (baseband) and `captures/run_*.csv` (log).
`run_identification.py` starts the horn counter-rotation and the capture
together; pass extra capture flags after `--`, e.g. `-- --gain 40 --ppm 1`.

## 4. Analyze

```
python3 GNURadioCode/analyzeall.py --fs 8000 --src-off 0 --pilot-off 3000 \
    0:captures/run_static.c64 2:captures/run_2rpm.c64 \
    4:captures/run_4rpm.c64 10:captures/run_10rpm.c64
```

This prints, per run, theory vs measured shift and the recovered radius, then
across runs the R^2, the linearity ratios, and whether the null is clean, and
writes the figures + `summary.csv`. For one run in detail use
`analyzecapture.py`.

## 5. Test the analysis without hardware

```
cd GNURadioCode
for r in 0 2 4 10; do python3 genv4cap.py --rpm $r --fs 8000 --src-off 0 --pilot-off 3000 --out m_${r}rpm; done
python3 analyzeall.py --fs 8000 --src-off 0 --pilot-off 3000 0:m_0rpm.c64 2:m_2rpm.c64 4:m_4rpm.c64 10:m_10rpm.c64
```

A matched example set is already under `GNURadioCode/simtest/`.

## What I need from you to finish the hardware wiring

The software path is complete and the analysis is verified on simulated data.
These values are placeholders in the code and must match your build before a
real run - tell me each and I will set them (and fix anything that is wrong):

1. **Central motor interface.** `dopplerwaterfall.py` assumes a Modbus-RTU BLDC
   driver on `/dev/ttyUSB0` at registers `0x8106` (control), `0x8110` (speed),
   `0x8206` (actual). Confirm your driver (BLD-515C?) uses these registers and
   that the speed register takes rpm directly. If it is PWM/analog instead, that
   path needs replacing.
2. **Horn motor wiring.** `horn_sync.py` / `MotorTools.py` default to BCM pins
   fwd 17, back 18, encoder A 22, B 23, and `counts_per_rev = 2048`. Give me the
   real H-bridge pins, encoder pins, and the goBILDA encoder counts per output
   rev. The PID gains will also need bench tuning.
3. **ADF4351 wiring.** SPI bus/dev (default 0/0), LE on GPIO 25, lock-detect on
   GPIO 24. Confirm these and that MUXOUT is wired to the lock pin. Channel
   spacing is 10 kHz, which puts the pilot on a ~5 kHz grid; if you need it at a
   precise 3 kHz offset, say so and I will lower the spacing.
4. **Capture config for your real runs.** Baseband rate (code default
   `fs_bb = 8000`; your old sim files were 12000 - pick one), source offset
   (0?), pilot offset (3000?), SDR raw rate (`fs_raw`, default 1.024e6; the paper
   used 2.048e6), center frequency, gain, ppm.
5. **Pointing.** Is manual aim-at-the-Sun fine for these runs (the paper's were
   ~10 min at the Sun), or do you want `Main.py` wired to the motors for real
   auto-tracking? `applyMotorStep` is currently a stub.
6. **The committed `.c64/.csv`.** They are simulated and inconsistent. OK to
   treat as throwaway? I left a clean matched set in `simtest/`.
