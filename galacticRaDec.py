# Find the RA and Dec crossing points of the galactic plane
# HISTORY
# 24Oct15 GIL parameterize number of samples 
# 24Feb28 GIL add checks and write output 
# 24Feb27 GIL major re=write;
# 21JAN15 GIL initial version based on web:
# https://learn.astropy.org/rst-tutorials/Coordinates-Transform.html
# Third-party dependencies

import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

# Set up matplotlib and use a nicer set of plot parameters
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
plt.style.use(astropy_mpl_style)
fig, ax = plt.subplots()

nargs = len( sys.argv)
#no arguments, really a no-op
if nargs < 2:
    print("galacticRaDec: computes Ra, and Dec for galacitic plane corrosing points{")
    print("usage")
    print("python galacticRaDec.py <outfilename>")
    print("where:")
    print("<outfilename> is the name of the galactic crossing file")
    print("")
    print("Glen Langston NSF - 24 February 27")
    exit()
outname = sys.argv[1]

ras = []
decs = []
# set the maximum number of galactic RA/Dec samples to compute
NSAMPLES = 1000
lons = np.linspace( 0., 360., NSAMPLES)
lat = 0.0
          
for lon in lons:
    c = SkyCoord(lon*u.degree, lat*u.degree, frame='galactic')
    ras.append( c.icrs.ra.deg)
    decs.append( c.icrs.dec.deg)
#    print(" %7.2f %7.2f" % (c.icrs.ra.deg,c.icrs.dec.deg))

ras = np.asarray(ras)
decs = np.asarray(decs)

ndec = len(decs)
declist = np.linspace(-90., 90, NSAMPLES)
nlist = len(declist)
# keep closest crossing ras and longitudes
ra1list = np.linspace(0., 360., NSAMPLES)
ra2list = np.linspace(0., 360., NSAMPLES)
lon1list = np.linspace(0., 360., NSAMPLES)
lon2list = np.linspace(0., 360., NSAMPLES)

#  get very accurate coordinates of galactic north and south poles
spole = SkyCoord(0.*u.deg, -90.*u.deg, frame='galactic')  # using degrees directly
npole = SkyCoord(0.*u.deg, +90.*u.deg, frame='galactic')  # using degrees directly

nPole = npole.icrs.ra.deg
nPoleDec = npole.icrs.dec.deg
sPole = spole.icrs.ra.deg
sPoleDec = spole.icrs.dec.deg

print("RA, Dec of North Galactic Pole: %9.5f, %9.5f" % \
      (nPole, npole.icrs.dec.deg))
print("RA, Dec of South Galactic Pole: %9.5f, %9.5f" % \
      (sPole, spole.icrs.dec.deg))

# for all decs in range -90 to 90
# find the closest 
id1 = 0
for iii in range(nlist):
    dec = declist[iii]
    mindec1 = 180.
    # find galactic latitude, longitude with minimum offset
    for ii in range( ndec):
        # the RAs great than north pole are in the 2nd list
        if not (ras[ii] < nPole and ras[ii] > sPole):
            continue
        ddec = np.abs(dec - decs[ii])
        if ddec < mindec1:
            id1 = ii
            mindec1 = ddec
    ra1list[iii] = ras[id1]
    # all lat 1s are zero
    lon1list[iii] = lons[id1]
    
# now should have two RAs for every declination
id2 = ndec-1
for iii in range(nlist):
    dec = declist[iii]
    mindec2 = 180.
    # find galactic latitude, longitude with minimum offset
    for ii in range( ndec):
        # the RAs great than north pole are in the 2nd list
        ddec = np.abs(dec - decs[ii])
        if ras[ii] < nPole and ras[ii] > sPole:
            continue
        if ddec < mindec2:
            id2 = ii
            mindec2 = ddec
    ra2list[iii] = ras[id2]
    lon2list[iii] = lons[id2]

f = open( outname, "w")
f.write( "# Table of Angular offsets between Transit Crossings\n")
f.write( "# of the Galactic Plane\n")
f.write( "# Glen Langston -- NSF -- glangsto@nsf.gov\n")
f.write( "#  Dec     RA1     RA2   Delta RA  Lon1    Lon2\n")
f.write( "#  (d)     (d)     (d)     (d)      (d)     (d) \n")

print("#N     Dec      Ra1     Ra2     Lon1    Lon2       ")
print("#      (d)      (d)     (d)      (d)     (d)       ")
# now save results and show a few values
for iii in range(nlist):
    # every so often print
    if iii == int(iii/30)*30:
        print("%3d %7.2f: %7.2f %7.2f %7.1f %7.1f" % \
          (iii, declist[iii], ra1list[iii], ra2list[iii], \
           lon1list[iii], lon2list[iii]))
    f.write("%7.2f %7.2f %7.2f %7.2f %7.1f %7.1f \n" % \
          (declist[iii], ra1list[iii], ra2list[iii], \
           ra2list[iii] - ra1list[iii], lon1list[iii], lon2list[iii]))
    

f.close()

plt.plot( ras, decs, label="RA,Dec", linewidth=5)
plt.plot( ra1list, declist, label="1st Crossing", color="cyan")
plt.plot( ra2list, declist, label="2nd Crossing")
plt.plot( ra2list - ra1list, declist, label="dRa", linestyle="dashed")
plt.plot([nPole],[nPoleDec], marker="+")
plt.plot([sPole],[sPoleDec], marker="x")

plt.annotate("  Galactic North", [nPole,nPoleDec])
plt.annotate("%7.3f,%7.3f " % (nPole, nPoleDec), [nPole,nPoleDec-10.])
plt.annotate("%7.3f,%7.3f " % (sPole, sPoleDec), [sPole,sPoleDec-10.])

plt.annotate("  Galactic South", [sPole,sPoleDec])
plt.xlabel("RA (deg)")
plt.ylabel("Dec (deg)")
ax.legend()
plt.title('Ra and Dec of Galactic Latitude crossing Points')
plt.tight_layout()
plt.show()
