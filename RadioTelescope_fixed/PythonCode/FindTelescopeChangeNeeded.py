import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.time import Time


# ---------- Astronomical helper functions ----------

def local_sidereal_time(jd, longitude_deg):
    JD0 = 2451545.0
    D = jd - JD0
    GMST = 280.46061837 + 360.98564736629 * D
    GMST = GMST % 360
    LST = (GMST + longitude_deg) % 360
    return LST


def hour_angle(lst_deg, ra_hours):
    return (lst_deg - ra_hours * 15) % 360


def radec_to_altaz(ra_hours, dec_deg, lat_deg, lon_deg, dt):
    jd = Time(dt, scale='utc').jd
    lst = local_sidereal_time(jd, lon_deg)
    ha = hour_angle(lst, ra_hours)

    ha = np.deg2rad(ha)
    dec = np.deg2rad(dec_deg)
    lat = np.deg2rad(lat_deg)

    alt = np.arcsin(
        np.sin(dec) * np.sin(lat) +
        np.cos(dec) * np.cos(lat) * np.cos(ha)
    )

    az = np.arctan2(
        -np.sin(ha),
        np.tan(dec) * np.cos(lat) - np.sin(lat) * np.cos(ha)
    )

    alt = np.rad2deg(alt)
    az = (np.rad2deg(az) + 360) % 360

    return alt, az


# ---------- Telescope Class ----------

class TelescopeOperations:

    def __init__(self, currLatitude, currLongitude, currElevation_km, desiredObject, currDateTime):
        self.currLatitude = currLatitude
        self.currLongitude = currLongitude
        self.currElevation = currElevation_km
        self.desiredObject = desiredObject

        # Telescope state (mount coordinates)
        self.currAlt = 90.0   # start pointing straight up
        self.currAz = 0.0     # facing north

        self.desiredAlt = 90.0
        self.desiredAz = 0.0

        self.updateDesiredValues(currDateTime)

    # In FindTelescopeChangeNeeded.py

    def updateDesiredValues(self, currDateTime):
        # Convert datetime to Julian Date for the API
        jd_time = Time(currDateTime).jd

        obj = Horizons(
            id=self.desiredObject,
            location={
                'lon': self.currLongitude,
                'lat': self.currLatitude,
                'elevation': self.currElevation
            },
            epochs=jd_time
        )

        ephem = obj.ephemerides()
        ra = ephem['RA'][0]
        dec = ephem['DEC'][0]

        alt, az = radec_to_altaz(
            ra, dec,
            self.currLatitude,
            self.currLongitude,
            currDateTime
        )

        self.desiredAlt = alt
        self.desiredAz = az

    def findAltAzChange(self):
        dAlt = self.desiredAlt - self.currAlt
        dAz = (self.desiredAz - self.currAz + 180) % 360 - 180
        return dAlt, dAz

    def applyMotorStep(self, dAlt, dAz, gain=0.25):
        self.currAlt += gain * dAlt
        self.currAz += gain * dAz
