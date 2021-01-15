# Find the Current Az and El of a few galactic coordinates
# HISTORY
# 21JAN15 GIL initial version based on web:
# https://learn.astropy.org/rst-tutorials/Coordinates-Transform.html
# Third-party dependencies
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
# next get components needed to compute Az,El
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.coordinates import AltAz

# Set up matplotlib and use a nicer set of plot parameters
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
plt.style.use(astropy_mpl_style)
fig, ax = plt.subplots()

# %matplotlib inline

#hcg7_center = SkyCoord(9.81625*u.deg, 0.88806*u.deg, frame='icrs')  # using degrees directly
#print(hcg7_center)

nargs = len( sys.argv)
if nargs < 2:
    print("toazel: computes Elevation of selected galactic coordinates for the input location")
    print("usage")
    print("toazel [<location>]")
    print("where <location> is the string identifying your telescope location")
    print("ie:  toazel GBT")
    print("or")
    print("toazel longitude latitude")
    print("where longitude is the telescope longitude in degrees (negative for west)")
    print("      latitude is the telescope latitude in degrees")
    print("Know locations are:")
    print(EarthLocation.get_site_names())
    print("")
    print("Glen Langston NSF - 21 January 15")
    exit()

if nargs == 2:
    mysite = sys.argv[1]
    print("Computing Elevations for site %s" % (mysite))
# Read GBT location; very close to my location.  For your location use the above
    asite = EarthLocation.of_site(mysite)
    print(asite)
else:
    try:
        lon=float( sys.argv[1])
    except:
        print( "Can Not parse argument %s as a floating point telescope longitude" % (sys.argv[1]))
        exit()
    try:
        lat=float( sys.argv[2])
    except:
        print( "Can Not parse argument %s as a floating point telescope latitude" % (sys.argv[2]))
        exit()
    asite = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=1000*u.m)
    print(asite)
    

tangent1 = SkyCoord(45.*u.deg, 0.*u.deg, frame='galactic')  # using degrees directly
tangent2 = SkyCoord(-45.*u.deg, 0.*u.deg, frame='galactic')  # using degrees directlyn
center   = SkyCoord(0.*u.deg, 0.*u.deg, frame='galactic')  # using degrees directly
anticenter = SkyCoord(180.*u.deg, 0.*u.deg, frame='galactic')  # using degrees directly
npole = SkyCoord(0.*u.deg, 90.*u.deg, frame='galactic')  # using degrees directly
spole = SkyCoord(0.*u.deg, -90.*u.deg, frame='galactic')  # using degrees directly

print(tangent1)
print(tangent1.icrs)

# Kitt Peak, Arizona
kitt_peak = EarthLocation(lat='31d57.5m', lon='-111d35.8m', height=2096*u.m)

#observing_time = Time('2010-12-21 1:00')
#aa = AltAz(location=kitt_peak, obstime=observing_time)

#print(aa)

#hcg7_center.transform_to(aa)

now = Time.now()

aa = AltAz(location=asite, obstime=now)
tangent1.transform_to(aa)

# this gives a Time object with an *array* of times
delta_hours = np.linspace(0, 24, 600)*u.hour
full_night_times = now + delta_hours
full_night_aa_frames = AltAz(location=asite, obstime=full_night_times)
full_night_aa_coos = tangent1.transform_to(full_night_aa_frames)


plt.plot(delta_hours, full_night_aa_coos.alt, label="l,b=45,0")

imax = np.argmax( full_night_aa_coos.alt)
az = full_night_aa_coos.az[imax].value
el = full_night_aa_coos.alt[imax].value + 1
hour = delta_hours[imax].value
t = ax.text( hour, el, "Az=%4.0f" % az)
#print( "Max Az %5.0f at index %d, hour=%s %5.0f " % ( az, imax, hour, el ))

full_night_aa_coos = tangent2.transform_to(full_night_aa_frames)

plt.plot(delta_hours, full_night_aa_coos.alt, label="l,b=-45,0")
imax = np.argmax( full_night_aa_coos.alt)
az = full_night_aa_coos.az[imax].value
el = full_night_aa_coos.alt[imax].value + 1
hour = delta_hours[imax].value
t = ax.text( hour, el, "Az=%4.0f" % az)

full_night_aa_coos = center.transform_to(full_night_aa_frames)

plt.plot(delta_hours, full_night_aa_coos.alt, linestyle='--', label="l,b=0,0")
imax = np.argmax( full_night_aa_coos.alt)
az = full_night_aa_coos.az[imax].value
el = full_night_aa_coos.alt[imax].value + 1
hour = delta_hours[imax].value
t = ax.text( hour, el, "Az=%4.0f" % az)

full_night_aa_coos = anticenter.transform_to(full_night_aa_frames)

plt.plot(delta_hours, full_night_aa_coos.alt, linestyle='-.', label="l,b=180,0")
imax = np.argmax( full_night_aa_coos.alt)
az = full_night_aa_coos.az[imax].value
el = full_night_aa_coos.alt[imax].value + 1
hour = delta_hours[imax].value
t = ax.text( hour, el, "Az=%4.0f" % az)

full_night_aa_coos = npole.transform_to(full_night_aa_frames)

plt.plot(delta_hours, full_night_aa_coos.alt, linestyle='dotted', label="l,b=0,90")
imax = np.argmax( full_night_aa_coos.alt)
az = full_night_aa_coos.az[imax].value
el = full_night_aa_coos.alt[imax].value + 1
hour = delta_hours[imax].value
t = ax.text( hour, el, "Az=%4.0f" % az)

full_night_aa_coos = spole.transform_to(full_night_aa_frames)

plt.plot(delta_hours, full_night_aa_coos.alt, linestyle='dotted', label="l,b=0,-90")
imax = np.argmax( full_night_aa_coos.alt)
az = full_night_aa_coos.az[imax].value
el = full_night_aa_coos.alt[imax].value + 1
hour = delta_hours[imax].value
t = ax.text( hour, el, "Az=%4.0f" % az)

nowstr = str( now)
nowparts = nowstr.split('.')
nowstr = nowparts[0]
plt.xlabel(('Hours from now (%s UTC)' % nowstr))
plt.ylabel('Elevation (deg)')
plt.legend(loc="upper right")
#plt.ylim(0.9,3)
plt.xlim(0,24)
# prepare to print location of site for calculations
geo = asite.geodetic
#print( geo.lon)
lon = str(geo.lon)
# trim out decimal point to simplify label
parts = lon.split('.')
lon = parts[0]
lat = str(geo.lat)
parts = lat.split('.')
lat = parts[0]
if nargs == 2:
    sitestr = "%s lon,lat=%s,%s" % (mysite, lon, lat)
else:
    sitestr = "Site lon,lat=%s,%s" % (lon, lat)

plt.title('Path of Galactic Coordinates in Elevation \n %s' % (sitestr))
plt.tight_layout()
plt.show()

