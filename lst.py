"""
lst is used to calculate the Local Sideral Time based on
the input location

Glen Langston,  2016 August 21
"""

##################################################
# Imports
##################################################
import datetime
import angles
import ephem

def calclst(longitude=-79.8397, latitude=34.4331, elevation=0.):
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
    print('LST = %s' % (lst))
    return lst

calclst()
