from time import sleep
from datetime import datetime, timezone

from FindTelescopeChangeNeeded import TelescopeOperations


# ---------- User Input ----------

currLongitude = float(input("Current Longitude (deg, +E): "))
currLatitude = float(input("Current Latitude (deg, +N): "))
currElevation = float(input("Current Elevation (km): "))
desiredObject = input("Desired object ID (e.g. Moon, Mars, 399): ")

# ---------- Initialize ----------

currDateTime = datetime.now(timezone.utc)

telescope = TelescopeOperations(
    currLatitude,
    currLongitude,
    currElevation,
    desiredObject,
    currDateTime
)

# ---------- Tracking Loop ----------

print("\nTracking started (Ctrl+C to stop)\n")

try:
    while True:
        currDateTime = datetime.now(timezone.utc)
        telescope.updateDesiredValues(currDateTime)

        dAlt, dAz = telescope.findAltAzChange()
        telescope.applyMotorStep(dAlt, dAz)

        print(
            f"Target Alt/Az: {telescope.desiredAlt:7.3f}°, {telescope.desiredAz:7.3f}° | "
            f"Scope Alt/Az: {telescope.currAlt:7.3f}°, {telescope.currAz:7.3f}° | "
            f"ΔAlt: {dAlt:6.3f}°, ΔAz: {dAz:6.3f}°"
        )

        sleep(1)

except KeyboardInterrupt:
    print("\nTracking stopped.")
