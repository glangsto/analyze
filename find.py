#Python Script to read in all spectra and compute on-off spectra for a range
#HISTORY
#17SEP27 GIL work in RA-Dec coordinates
#17SEP26 GIL work in RA-Dec coordinates
#17SEP25 GIL init version based on m.py
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import statistics
import angles
import ephem
#import radec2azel
import radioastronomy
import copy
#from scipy.signal import savgol_filter
import interpolate

linelist = [1420.0, 1419.0, 1418.0]  # RFI lines in MHz
linewidth = [9, 5, 7]
#for 2017-07-21 observations
eloffset = -4.91
azoffset = -18.0
#for 2017-09-03 observations
eloffset = 0.52
azoffset = -15.36
fwhm = 16.0 # telescope fwhm angle degrees
onTheta = angles.d2r(0.4*fwhm)
offThetaA = angles.d2r(0.5*fwhm)
offThetaB = angles.d2r(1.0*fwhm)

nargs = len(sys.argv)

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1e8
xallmax = -9.e9
xallmin = 9.e9
yallmax = -9.e9
yallmin = 9.e9
# velocities for fitting baselines
minvel = -450.
minvel = -280.
maxvel = 210.
# currently used
maxvel = 160.
minvel = -160.
# currently used
maxvel = 180.
minvel = -180.

# define the reference location and FWHM
raStr = sys.argv[1]
decStr = sys.argv[2]

# will use sum within FWHM and take the reference as 1 to 2 FWHM
radec2000 = ephem.Equatorial( raStr, decStr, epoch=ephem.J2000)
print 'Ra,Dec = %s, %s (J2000)' % (radec2000.ra, radec2000.dec)
#radec = ephem.Equatorial(radec2000.ra, radec2000.dec, epoch=datestr)
gal = ephem.Galactic(radec2000)
print 'Lon,Lat= %s, %s (Galactic)' % (gal.lon, gal.lat)
# need to turn dd:mm:ss into float
aparts = angles.phmsdms(raStr)
ra0 = angles.sexa2deci(aparts['sign'], *aparts['vals'], todeg=True)
aparts = angles.phmsdms(decStr)
dec0 = angles.sexa2deci(aparts['sign'], *aparts['vals'], todeg=False)
print 'ra,dec: ', ra0,dec0
# convert to radians
ra0 = angles.d2r( ra0)
dec0 = angles.d2r( dec0)

c = 299792.  # (v km/sec)
nuh1 = 1420.4056 # neutral hydrogen frequency (MHz)
thot = 285.0  # define hot and cold
#thot = 272.0  # 30 Farenheit = 272 K
tcold = 10.0
tmin = 20.0 
tmax = 999.0 # define reasoanable value limits

# prepare to remove a linear baseline

xa = 200
#xb = 1024-xa
xb = 1024-200

#for symbol, value in locals().items():
#    print symbol, value

nplot = 0
nhot = 0         # number of obs with el < 0
ncold = 0
nhigh = 0        # highest galactic latitude
minGlat = +90.
maxGlat = -90.
lastfreq = 0.
lastbw = 0.
lastgain = 0.
lastel = 0.
lastaz = 0.
firstdate = ""
lastdate = ""
minel = 200.
maxel = -200.

# simple argument parsing
doFlag = True
firstfile = 3
# initialize the minimum distance (radians) to a big value
minTheta = 3.1415926
iTheta = 0

# first read through all data and find hot load
names = sys.argv[firstfile:]
names = sorted(names)
nfiles = len(names)
nData = 0
nRead = 0
nOn = 0
nOff = 0

for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
# only working with positive elevations
    if rs.telel < 0:
        continue
    rs.telel = rs.telel + eloffset
    rs.telaz = rs.telaz + azoffset
    rs.azel2radec()    # compute ra,dec from az,el
    
    # if spectra array is not yet initialized
    if nData == 0:
        nData = len(rs.xdata)
        spectra = np.zeros((nfiles,nData))
        ras = np.zeros((nfiles))
        decs = np.zeros((nfiles))
        azs = np.zeros((nfiles))
        els = np.zeros((nfiles))
        thetas = np.zeros((nfiles))
        utcs = np.zeros((nfiles))
        rsmin = copy.deepcopy( rs)
        rsmax = copy.deepcopy( rs)
        rsmed = copy.deepcopy( rs)
        rsave = copy.deepcopy( rs)
        rsave = copy.deepcopy( rs)
        xv = rs.xdata*1.E-6
        print "Telescope Az, El: ", rs.telaz, rs.telel

# convert to MHz
    yv = rs.ydataA * scalefactor
    if doFlag:
        hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
        spectra[nRead, : ] = hv
    else:
        spectra[nRead, : ] = yv
    # now record the coordinates to compare distances.
    azs[nRead] = rs.telaz
    els[nRead] = rs.telel
    ras[nRead] = rs.ra
    decs[nRead] = rs.dec
#    print rs.ra, rs.dec
    aTheta = angles.sep(angles.d2r(rs.ra), angles.d2r(rs.dec), ra0, dec0)
    thetas[nRead] = aTheta
    # keep track of minimum angular distance
    if aTheta < minTheta:
        minTheta = aTheta
        iTheta = nRead

    if aTheta < onTheta:
        if nOn == 0:
            onSpectra = copy.deepcopy( rs)
        else:
            onSpectra.ydataA = onSpectra.ydataA + rs.ydataA
        nOn = nOn + 1
        nRead = nRead + 1

    if aTheta > offThetaA and aTheta < offThetaB:
        if nOff == 0:
            offSpectra = copy.deepcopy( rs)
        else:
            offSpectra.ydataA = offSpectra.ydataA + rs.ydataA
        nOff = nOff + 1
        nRead = nRead + 1

#    utcs[nRead] = rs.utc
#    lsts[nRead] = rs.lst


print 'Min Theta: ', angles.r2d( minTheta)
print 'N On: ', nOn, ' N Off: ', nOff

print 'Read ',nRead,' Spectra with ',nData,' Spectral channels '

if nOn < 1:
    print 'no On Spectra found'
    exit()

onSpectra.ydataA = scalefactor*onSpectra.ydataA/float(nOn)

if nOff < 1:
    print 'no Off Spectra found'
    exit()

offSpectra.ydataA = scalefactor*offSpectra.ydataA/float(nOff)

deltaSpectra = copy.deepcopy(onSpectra)
deltaSpectra.ydataA = deltaSpectra.ydataA-offSpectra.ydataA

for iii in range( 0, nData):
    if offSpectra.ydataA[iii] > 0.1:
        deltaSpectra.ydataA[iii] = 100.*deltaSpectra.ydataA[iii]/offSpectra.ydataA[iii]
    else:
        deltaSpectra.ydataA[iii] = 0

fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)

#plt.plot(xv, onSpectra.ydataA, colors[0], linestyle=linestyles[0],label='ON', lw=4)
#plt.plot(xv, offSpectra.ydataA, colors[1], linestyle=linestyles[1],label='OFF', lw=4)
plt.plot(xv, deltaSpectra.ydataA, colors[1], linestyle=linestyles[1],label='Delta', lw=4)

# now do channel by channel statistics
for iii in range( 0, nData):
    ydata = spectra[:,iii]
    ymed = statistics.median(ydata)
    ymin = min(ydata)
    ymax = max(ydata)
    yave = statistics.mean(ydata)
    rsmin.ydataA[iii] = ymin
    rsmax.ydataA[iii] = ymax
    rsmed.ydataA[iii] = ymed
    rsave.ydataA[iii] = yave

# now find the min and max spectra by name
# compare values near the middle but not exacly at the middle
middle = int(nData/2) - 20 
imin = 0
imax = 0
amin = spectra[imin,middle]
amax = spectra[imin,middle]

# now for all other files compare intensities of middle values
for iii in range( 1, nRead):
    if spectra[imin,middle] > spectra[iii,middle]:
        imin = iii
    if spectra[imax,middle] < spectra[iii,middle]:
        imax = iii
    
#fig, ax1 = plt.subplots(figsize=(10, 6))

#plt.plot(xv, rsmax.ydataA, colors[1], linestyle=linestyles[1],label='Maximum', lw=4)
#plt.plot(xv, rsmed.ydataA, colors[2], linestyle=linestyles[2],label='Median', lw=4)
#plt.plot(xv, rsmin.ydataA, colors[0], linestyle=linestyles[0],label='Minimum', lw=4)
#plt.plot(xv, rsave.ydataA, colors[3], linestyle=linestyles[3],label='Average', lw=4)

print 'Maximum Spectrum File:'
print names[imax]

rs.read_spec_ast(names[imax])
gallon = rs.gallon
gallat = rs.gallat
label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( names[imax],rs.telaz,rs.telel,gallon,gallat)
#plt.plot(xv, rs.ydataA*scalefactor, colors[5], linestyle=linestyles[5],label=label, lw=4)
rs.read_spec_ast(names[imin])
gallon = rs.gallon
gallat = rs.gallat
label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( names[imin],rs.telaz,rs.telel,gallon,gallat)
#plt.plot(xv, rs.ydataA*scalefactor, colors[4], linestyle=linestyles[4],label=label, lw=4)
print 'Minimum Spectrum File:'
print names[imin]

mytitle = '%s - %s' % ( names[0], names[nRead-1])
fig.canvas.set_window_title(mytitle)

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

plt.title(mytitle, fontsize=16)
plt.xlabel('Frequency (MHz)', fontsize=16)
plt.ylabel('Intensity (Counts)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
