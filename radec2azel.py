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

def calcazel(ra, dec, longitude=-79.8397, latitude=38.4331, elevation=800.):
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
    fmt = 'Date   = %s,  LST = %s'
    print fmt % (datestr, lst)
    # first put ra,dec in ephemeris format
    radec2000 = ephem.Equatorial( ra, dec, epoch=ephem.J2000)
    # now compute ra, dec of date
    radec = ephem.Equatorial(radec2000.ra, radec2000.dec, epoch=datestr)
    print 'Ra,Dec = %s, %s (J2000)' % (radec2000.ra, radec2000.dec)
    gal = ephem.Galactic(radec2000)
    print 'Lon,Lat= %s, %s (Galactic)' % (gal.lon, gal.lat)
    object = ephem.FixedBody()
    object._ra = radec.ra
    object._dec = radec.dec
    object.compute(location)
    # need to turn dd:mm:ss into float
    aparts = angles.phmsdms(str(object.az))
    az = angles.sexa2deci(aparts['sign'], *aparts['vals'], todeg=False)
    aparts = angles.phmsdms(str(object.alt))
    el = angles.sexa2deci(aparts['sign'], *aparts['vals'], todeg=False)
#    print "Az: %s -> %f" % (str(object.az), az)
#    print "El: %s -> %f" % (str(object.alt), el)
    fmt = 'Az, El = %.2f, %.2f (degrees)' 
    print fmt % (az, el)
    
    return az, el

#print 'Argument List:', str(sys.argv)
nargs = len( sys.argv)

if nargs < 3:
    print 'radec2azel: compute the Ra,Dec (2000) from input az,el (degrees)'
    exit()

calcazel( sys.argv[1], sys.argv[2])

