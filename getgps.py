#! /usr/bin/python
#Python Script to average GPS values and return estimate of uncertainty.
#HISTORY
# 21Jul08 GIL print arguments for FIX 
# 21Apr08 GIL return rms values 
# 21Apr02 GIL initial version

import sys
from gps import *
import time
import numpy as np

MAXCOUNT = 10
maxcount = MAXCOUNT
nargs = len( sys.argv)
# assume first argument is number of samples to average:
if nargs > 1:
    try:
        maxcount = int(sys.argv[1])
    except:
        print("Invalid number of samples to average: %s" % (sys.argv[1]))


print("Averaging %d GPS measurements" % (maxcount))
tellons = np.zeros( maxcount)
tellats = np.zeros( maxcount)
telalts = np.zeros( maxcount)

count = 0
gpsd = gps(mode=WATCH_ENABLE|WATCH_NEWSTYLE) 
print('latitude\tlongitude\ttime utc\t\t\taltitude\tepv\tept\tspeed\tclimb')
# '\t' = TAB to try and output the data in columns.
avelat = 0.
avelon = 0.
avealt = 0.

try:
 
    while count < maxcount:
        report = gpsd.next() #
        if report['class'] == 'TPV':
             
            print ( "%s %s" % (getattr(report,'lat',0.0),"\t"), end="")
            print ( "%s %s" % (getattr(report,'lon',0.0),"\t"), end="")
            print ( "%s %s" % (getattr(report,'time',''),"\t"), end="")
            print ( "%s %s" % (getattr(report,'alt','nan'),"\t"), end="")
            print ( "%s %s" % (getattr(report,'epv','nan'),"\t"), end="")
            print ( "%s %s" % (getattr(report,'ept','nan'),"\t"), end="")
            print ( "%s %s" % (getattr(report,'speed','nan'),"\t"), end="")
            print( "%s %s" % (getattr(report,'climb','nan'),"\t")
)
            # now sum averages
            avelon = avelon + float(getattr(report,'lon',0.0))
            tellons[count] = float(getattr(report,'lon',0.0))
            avelat = avelat + float(getattr(report,'lat',0.0))
            tellats[count] = float(getattr(report,'lat',0.0))
            avealt = avealt + float(getattr(report,'alt',0.0))
            telalts[count] = float(getattr(report,'alt',0.0))
            count = count + 1
        time.sleep(1) 
    avelon = avelon / float(count)
    avelat = avelat / float(count)
    avealt = avealt / float(count)
    dlon = np.std( tellons)
    dlat = np.std( tellats)
    dalt = np.std( telalts)
    print("Average lat, lon, alt: %13.9f %13.9f %9.3f" % \
          ( avelat, avelon, avealt))
    print("   +/-  lat, lon, alt: %13.9f %13.9f %9.3f" % \
          ( dlat, dlon, dalt))
    rEarth = 6.371E6
    dlat = dlat * rEarth * np.pi / 180.
    dlon = dlon * rEarth * np.pi / 180.
    print("   +/-m lat, lon, alt: %10.3f %10.3f %9.3f (meters)" % \
          ( dlat, dlon, dalt))
    print("FIX arguments for updating notes files:")
    print("FIX -RE -LAT %11.7f -LON %11.7f -ALT %7.3f" % \
          (avelat, avelon, avealt))

# keep latitude and longitude in float format
    lat = avelat
    lon = avelon

   # convert to degrees minutes seconds
   # get sign of latitude
    if lat >= 0:
        pmlat = '+'
    else:
        pmlat = '-'
        lat = - lat

    # get sign of longitude
    if lon >= 0:
        pmlon = '+'
    else:
        pmlon = '-'
        lon = - lon

# get degrees part
    dlat = int(lat)
    dlon = int(lon)

# start on minutes part
    lat = (lat - dlat)*60.
    lon = (lon - dlon)*60.

    mlat = int(lat)
    mlon = int(lon)

# start on seconds part
    slat = (lat - mlat)*60.
    slon = (lon - mlon)*60.

# show results
    print(('Latitude : %10.6f' % (avelat)))
    print(('Latitude : %s%02d:%02d:%05.2f' % (pmlat,dlat,mlat,slat)))
    print(('Longitude: %10.6f' % (avelon)))
    print(('Longitude: %s%02d:%02d:%05.2f' % (pmlon,dlon,mlon,slon)))

except (KeyboardInterrupt, SystemExit): #when you press ctrl+c
    print( "Done.\nExiting.")

