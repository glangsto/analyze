"""
lst is used to calculate the Local Sideral Time based on
the input location

Glen Langston,  2016 August 21
"""

##################################################
# Imports
##################################################
import sys
import datetime
import angles
import ephem

def calclst(longitude=-79.8397, latitude=34.4331, elevation=810.):
    """
    Compute the ra,dec (J2000) from Az,El location and time
    """
    location = ephem.Observer()
    location.lon = angles.d2r(longitude)
    location.lat = angles.d2r(latitude)
    location.elevation = elevation
    now = datetime.datetime.utcnow()
    strnow = now.isoformat()
    dates = strnow.split('T')
    datestr = dates[0] + ' ' + dates[1]
    location.date = datestr
    lst = location.sidereal_time()
    print 'LST = %s' % lst
    return lst

# Calculate LST time: input parameters are 
# longitude (degrees)  East is positive, West is negative
# latidute (degrees) North is positive, South is negative
# elevation (m above mean sea level)
# default values are for Green Bank, WV

nargs = len( sys.argv)

if nargs > 1:
    long = float(sys.argv[1])
    print "Computing LST for Longitude %9.6f (d)" % (long)
if nargs > 2:
    lat = float(sys.argv[2])
    print "Computing LST for Latitude  %9.6f (d)" % (lat)
if nargs > 3:
    el = float(sys.argv[3])
    print "Computing LST for Elevation %9.6f (m)" % (el)

if nargs == 1:
    calclst()
elif nargs == 2:
    calclst(longitude=long)
elif nargs == 3:
    calclst(longitude=long, latitude=lat)
else:
    calclst(longitude=long, latitude=lat, elevation=el)
