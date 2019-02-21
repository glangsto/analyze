#Python Script to create an RFI list
#plot the raw data from the observation
#HISTORY
#17MAR23 GIL Initial version
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import statistics
import radioastronomy
import copy
import findRfi

from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

import interpolate

avetimesec = 3600.
dy = -1.
nlines = 0
MAXLINES = 100
lineList = np.zeros(MAXLINES)
lineChannel = np.zeros(MAXLINES)
lineWidth = np.zeros(MAXLINES)

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
minvel = -160.
maxvel = 160.

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

#first argument is the averaging time in seconds

avetimesec = float(sys.argv[1])
print "Average time: ", avetimesec, " (seconds)"

# first read through all data and find hot load
names = sys.argv[2:]
names = sorted(names)
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el

    if rs.telel < 0:
        if nhot == 0:
            hot = copy.deepcopy( rs)
            nhot = 1
        else:
            hot.ydataA = hot.ydataA + rs.ydataA
            hot.count = hot.count + rs.count
            nhot = nhot + 1
    else: # else above horizon, find min, max galactic latitudes
        if rs.gallat > maxGlat:
            maxGlat = rs.gallat
        if rs.gallat < minGlat:
            minGlat = rs.gallat

if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA / float(nhot)
    print "Found %3d hot load obs" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.E-6
yv = hot.ydataA
nData = len(xv)
#previously just copy
hv = np.zeros(nData)
for iii in range(nData):
    hv[iii] = yv[iii]

vel = np.zeros(nData)
for jjj in range (0, nData):
    vel[jjj] = c * (nuh1 - xv[jjj])/nuh1

print 'Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % (minGlat, maxGlat)

# assume only a limited range of galactic latitudes are available
# not range about +/-60.
use60Range = False

# all galactic latitudes above +/-60d can be used
if minGlat < -60. or maxGlat > 60.:
    minGlat = -60.
    maxGlat = 60.
else: # else no high galactic latitude data
    # use highest galactic latitudes - +/-5.degrees
    if -minGlat > maxGlat:  # if negative latitudes higher
        minGlat = minGlat + 5.
        maxGlat = 90.
    else: # else positive latitudes higher
        maxGlat = maxGlat - 5.
        minGlat = -90.
print "Min, Max Glat for cold load: ", minGlat, maxGlat

# now average coldest data for calibration
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el

    if rs.telel < 0:
        continue
    if rs.gallat > maxGlat or rs.gallat < minGlat:
        if nhigh == 0:
            high = copy.deepcopy( rs)
            high.ydataA = rs.ydataA
            nhigh = 1
        else:
            high.ydataA = high.ydataA + rs.ydataA
            high.count = high.count + rs.count
            nhigh = nhigh + 1

if nhigh < 1.:
    print "No high galactic latitude data: can not calibrate"
    exit()
else:
    high.ydataA = scalefactor * high.ydataA/nhigh
    print "Found %d High Galactic Latitude spectra" % (nhigh)
    yv = high.ydataA

cv = np.zeros(nData)
for iii in range(nData):
    cv[iii] = yv[iii]
#    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

nMedian = int(11)

hotValues, hotMedian = findRfi.findRfi( xv, hv, nMedian)

# prepare to take the median baseline
coldValues, coldMedian = findRfi.findRfi( xv, cv, nMedian)

# finally compute gain on a channel by chanel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hotMedian[iii] - coldMedian[iii])/(thot - tcold)

nEdge = 50
medianGain = statistics.median[nEdge:(nData-nEdge)]

trx = np.zeros(nData)
tsys = np.zeros(nData)
for iii in range(nData):
    trx[iii] = (coldMedian[iii]+coldValues[iii])*gain[iii]) - tcold

tsysMedian = statistics.median(trx)

for iii in range(nData):
    gain[iii] = (tsysMedian + thot)/hotMedian[iii]

fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax

plt.plot(xv, hotMedian, 'r', linestyle='-.',label="Hot")
plt.plot(xv, coldMedian, 'b', linestyle='-.',label="Cold")
plt.plot(xv, hotValues, 'r', linestyle='-',label="Hot Residual")
plt.plot(xv, coldValues, 'b', linestyle='-',label="Cold Residual")

hotRfi = hv - (hotMedian + hotValues) + 1.
coldRfi = cv - (coldMedian + coldValues) + 1.
plt.plot(xv, hotRfi, 'r', linestyle='--',label="Hot RFI Identified")
plt.plot(xv, coldRfi, 'b', linestyle='--',label="Cold RFI Identified")

plt.legend(loc='upper right')
plt.show()
# exit while median loop


Tsys = statistics.median(trx)
print "Tsys ", tsysMedian
plt.plot(xv, trx, 'b', linestyle='--',label="Trx")
plt.legend(loc='upper right')
plt.show()


