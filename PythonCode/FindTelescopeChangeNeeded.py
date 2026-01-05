# Take current location as if facing directly up and telescope oriented directly north
# From here, find desired right ascension and declination as well as hour angle and feed it in for precise locating
# If the body is far away, use function to find the additional movement needed.

class TelescopeOperations:
    currLatitude = 0
    currLongitude = 0
    desiredDeclination = 0
    desiredRightAscension = 0
    currDeclination = 0
    currRightAscension = 0
    def __init__(self, currLatitude, currLongitude, desiredDeclination, desiredRightAscension, currJulianDate):
        self.updateDesiredValues(desiredRightAscension,desiredDeclination)
        self.updateCurrentValues(currJulianDate, currLatitude, currLongitude)
    def updateCurrentValues(self,currJulianDate,currLatitude, currLongitude):
        self.currLatitude = currLatitude
        self.currLongitude = currLongitude
        self.currDeclination = self.currLatitude
        JD0 = 2451545.0
        D = currJulianDate - JD0
        thetaGMST = 280.46061837 + 360.98564736629 * D
        thetaGMST = divmod(thetaGMST, 360)
        thetaLST = thetaGMST + self.currLongitude
        if (thetaLST < 0):
            thetaLST += 360
        self.currRightAscension = thetaLST / 15
    def updateDesiredValues(self, desiredRightAscension, desiredDeclination):
        self.desiredRightAscension = desiredRightAscension
        self.desiredDeclination = desiredDeclination
    def findDeclinationChange(self):
        return self.desiredDeclination - self.currDeclination
    def findAscensionAngleChange(self):
        return self.desiredRightAscension * 15 - self.currRightAscension * 15
