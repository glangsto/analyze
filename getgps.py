#! /usr/bin/python
 
from gps import *
import time

MAXCOUNT = 10
count = 0
gpsd = gps(mode=WATCH_ENABLE|WATCH_NEWSTYLE) 
print 'latitude\tlongitude\ttime utc\t\t\taltitude\tepv\tept\tspeed\tclimb' # '\t' = TAB to try and output the data in columns.
avelat = 0.
avelon = 0.
avealt = 0.
try:
 
 
    while count < MAXCOUNT:
        report = gpsd.next() #
        if report['class'] == 'TPV':
             
            print  getattr(report,'lat',0.0),"\t",
            print  getattr(report,'lon',0.0),"\t",
            print getattr(report,'time',''),"\t",
            print  getattr(report,'alt','nan'),"\t\t",
            print  getattr(report,'epv','nan'),"\t",
            print  getattr(report,'ept','nan'),"\t",
            print  getattr(report,'speed','nan'),"\t",
            print getattr(report,'climb','nan'),"\t"
            # now sum averages
            avelon = avelon + float(getattr(report,'lon',0.0))
            avelat = avelat + float(getattr(report,'lat',0.0))
            avealt = avealt + float(getattr(report,'alt',0.0))
            count = count + 1
        time.sleep(1) 
    avelon = avelon / float(count)
    avelat = avelat / float(count)
    avealt = avealt / float(count)
    print("Average lat, lon, alt: %13.9f %13.9f, %9.3f" % \
          ( avelat, avelon, avealt))

except (KeyboardInterrupt, SystemExit): #when you press ctrl+c
    print "Done.\nExiting."
    x
