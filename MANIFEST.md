# Program manifest

Every code file in the project, what it does, and whether it is part of the
**identification** pipeline (capture a run -> recover the Doppler fingerprint).
Status key:

- **CORE** - needed for identification to work
- **UTILITY** - helper you will use but not strictly required
- **TEST** - bench / bring-up script, not used in a real run
- **LEGACY** - the older high-speed method, kept for reference, not identification
- **EMPTY** - file was lost / never filled in
- **NEW** - added now to close a gap

## GNURadioCode/  (pilot, capture, analysis)

| File | Status | In pipeline | What it does |
|------|--------|-------------|--------------|
| `pilot_adf4351.py` | CORE, NEW | yes | Pi -> ADF4351 pilot over SPI: register calc + program + lock check, plus a CLI. This is the Pi/ADF communication as its own module. |
| `dopplerwaterfall.py` | CORE | yes | The capture engine. Programs the pilot, drives the central motor (Modbus), records the decimated baseband `.c64`, logs `t/angle/rpm/src/pilot/corrected/SNR` to `.csv`, and tracks both tones live. Now imports `pilot_adf4351`. |
| `analyzecapture.py` | CORE | yes | Analyze one capture: drift-cancel, fit the sinusoid against the logged angle, print amplitude / velocity / recovered radius, save figures. |
| `analyzeall.py` | CORE | yes | Analyze a set of runs (0/2/4/10 rpm). Fits at the known rotation rate, writes waterfalls, cancellation and fit figures, the shift-vs-rpm chart, R^2, linearity, and `summary.csv`. Rewritten to take the config on the command line. |
| `genv4cap.py` | UTILITY | test aid | Simulate a capture (`.c64`+`.csv`) so the analysis can be tested without hardware. Rewritten with a CLI. |
| `waterfallmaker.py` | UTILITY | helper | Plot a spectrogram from any `.c64`. |
| `recordspectrum.py` | LEGACY | no | Old GNU Radio recorder: averaged FFT frames to `.dat`, stationary vs rotating. |
| `runthis.py` | LEGACY | no | Old peak shift / broadening comparison of two `.dat` files. |
| `SDRTest.py` | TEST | no | GNU Radio FM test flow graph (96.7 MHz). SDR bring-up. |
| `rtlsdrtest.py` | TEST | no | Real-time spectrum viewer at 1420 MHz. SDR bring-up. |
| `sdrtesting.py` | TEST | no | FM demod to audio. SDR learning / bring-up. |
| `main.py` | EMPTY | - | Was lost. The run entry point is now `run_identification.py` at the repo root. |
| `*.grc` | LEGACY/TEST | no | GNU Radio flow-graph sources for the test / legacy flows. |
| `*.c64 *.csv *.dat *.png` | DATA | - | Simulated / legacy data and figures. The committed `.c64/.csv` are simulated and **inconsistent** (some at source 0 / pilot 3000, some at 1500 / 4500). Treat as throwaway; a fresh matched set is in `simtest/`. |

## PythonCode/  (motors and pointing, on the Pi)

| File | Status | In pipeline | What it does |
|------|--------|-------------|--------------|
| `MotorTools.py` | CORE | yes | Horn-motor controller: encoder feedback + PID, angle interface. |
| `horn_sync.py` | CORE, NEW | yes | Counter-rotate the horn at -omega during a run so it stays on target. Standalone CLI or `run_horn_sync()` in a thread. |
| `FindTelescopeChangeNeeded.py` | SUPPORT | optional | Target Alt/Az from JPL Horizons; `TelescopeOperations`. Needed only for auto-tracking a non-Sun target. |
| `Main.py` | SUPPORT | optional | Pointing / tracking loop. Note: `applyMotorStep` is a stub that updates state but does not yet drive the motors. |
| `AccountForRotation.py` | SUPPORT | no | One-shot "move the horn by this angle". Superseded by `horn_sync.py` for a continuous run. |
| `RunSteppers.py` | TEST | no | Raw stepper pulse train. Central-motor bring-up. |
| `DCMotorRun.py` | TEST | no | L298N DC motor full-on. |
| `DCMotorTest.py` | TEST | no | L298N DC motor PWM ramp. |
| `StepperMotorTest.py` | TEST | no | Toggle a step pin to check wiring. |

## Repo root

| File | Status | In pipeline | What it does |
|------|--------|-------------|--------------|
| `run_identification.py` | CORE, NEW | yes | Master for one run: starts the horn counter-rotation and the capture together, stops both together. |
| `requirements.txt` / `apt-packages.txt` / `setup.sh` | SETUP | - | Dependencies and install. |

## The minimum identification pipeline

```
per run (on the Pi):   run_identification.py
                         |-- GNURadioCode/dopplerwaterfall.py   (pilot + central motor + SDR + log)
                         |      \-- GNURadioCode/pilot_adf4351.py
                         \-- PythonCode/horn_sync.py            (horn counter-rotation)
                                \-- PythonCode/MotorTools.py

offline:               GNURadioCode/analyzeall.py   (all runs -> figures + R^2 + summary.csv)
                       GNURadioCode/analyzecapture.py (one run in detail)

no hardware (testing): GNURadioCode/genv4cap.py -> analyzeall.py
```
