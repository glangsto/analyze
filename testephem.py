#Procedure to calculate the galatic latitude and logitude in an azimuth circle
#Assuming the telescope operates at fixed elevation.
#Glen Langston
#HISTORY
#16OCT10 GIL take input elevation (degrees)
#15JUL06 GIL add comments
#15JUN01 GIL initial version
import sys
from ephem import Observer
from ephem import Sun
from ephem import Equatorial
from ephem import Galactic
import datetime

# set the default azimuth and elevation for the  first calcuation
az = '220'
el = '20'
el = '45'

nargs = len(sys.argv)
el = float(sys.argv[1])

print 'Calculating galactic coordinates for elevation '+str(el)+' deg'

# first define the observer location
me = Observer()
me.lon='-79.8397'
me.lat='38.4331'
me.elevation=800   # height in meters
# reformat the time into ephemerus format. 
now = datetime.datetime.utcnow()
strnow = now.isoformat()
dates = strnow.split('T')
datestr = dates[0] + ' ' + dates[1]
me.date = datestr
timestrfull = str(dates[1])
times = timestrfull.split('.')
outname = dates[0] + '_' + times[0] + '.txt'
#print(outname)
# check the LST as a test
#print('UTC = ' + datestr + ' -> ' + str(me.sidereal_time()) + ' LST ')

#now compute the ra,dec from the azimuth,elevation, time and location
ra_a,dec_a = me.radec_of( az,str(el))
print('Horn Ra,Dec = %s,%s for Az,El %s,%s of %s' % (ra_a,dec_a,az,str(el), datestr))
#radec_a = Equatorial( ra_a, dec_a, me.date)
radec = Equatorial( ra_a, dec_a, epoch=datestr)
print(radec.ra, radec.dec)
#finally calculate the galactic coordiates for the ra, dec coordinates
gal = Galactic( radec)
print(gal.lon, gal.lat)

#now write a file with the coordinate sets
outfile = open( outname, 'w')
outline = '   Az       Ra          Dec          G-Lon        G-Lat \n'
outfile.write(outline)

#for all azimuths, assuming a fixed elevation
for azimuth in range(0, 360,10):
    az = str( azimuth) 
    ra_a,dec_a = me.radec_of( az, str(el))
    radec = Equatorial( ra_a, dec_a, epoch=datestr)
    gal = Galactic( radec)
    print( '%6s %12s %12s %12s %12s' % (az, ra_a, dec_a, gal.lon, gal.lat))
    outline = '%6s %12s %12s %12s %12s \n' % (az, ra_a, dec_a, gal.lon, gal.lat)
    outfile.write(outline)


# perform a sanity check by computeing the solar positions in az,el and lon,lat
sun = Sun(me)
outline = 'Sun = %12s,%12s (ra,dec)' % (sun.a_ra,sun.a_dec)
print(outline)
outfile.write(outline+'\n')

outline = 'Sun = %12s,%12s (alt,az)' % (sun.alt,sun.az)
print(outline)
outfile.write(outline+'\n')

sunradec = Equatorial( sun.a_ra, sun.a_dec, epoch=datestr)
sungal = Galactic( sunradec)

outline = 'Sun = %12s,%12s (lon,lat)' % (sungal.lon,sungal.lat)
print(outline)
outfile.write(outline+'\n')

outline = 'UTC = ' + datestr +' -> ' + str(me.sidereal_time()) + ' LST'
print(outline)
outfile.write(outline+'\n')

outfile.close()


