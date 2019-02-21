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
linelist = [1420.0, 1419.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 5, 7]

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
    print "Found %d High Galactic Latidue spectra" % (nhigh)
    yv = high.ydataA
#    yv[511] = ((3.*yv[510])+yv[514])/4.
#    yv[512] = (yv[510]+yv[514])/2.
#    yv[513] = ((yv[510])+(3.*yv[514]))/4.

#    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
    cv = np.zeros(nData)
    for iii in range(nData):
        cv[iii] = yv[iii]

# finally compute gain on a channel by chanel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)

fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
nMedian = 17

hotSignal, hotMedian = findRfi( xv, hv, nMedian)
# prepare to take the median baseline

nMedian = int(17)
myv = np.zeros(nData)
ys = np.zeros(nData)
xs = np.zeros(nData)
medianValues = np.zeros(nData)

nloop = 0
# loop through hot load removing RFI and writing RFI list
while nMedian > 7:

#    print 'N median: %d ' % (nMedian)
    for iii in range(nMedian,nData-(nMedian+1)):
        ys = hv[(iii-nMedian):(iii+nMedian)]
        myv[iii] = statistics.median(ys)
        medianValues[iii] = medianValues[iii] + myv[iii]

    dyv = np.zeros(nData)
    dyd = np.zeros(nData)
    dyv = hv - myv

    # flatten out the ends
    for iii in range(nMedian+1):
        dyv[iii] = 0
        dyv[nData-iii-1] = 0

    ymed = statistics.median(yv)
    yrms = statistics.stdev(yv)
    hmed = statistics.median(hv)
    hrms = statistics.stdev(hv)
    dyrms = statistics.stdev(dyv)
    count = hot.count
    ncold = 0

    print "RMS before and after median: ", hrms, dyrms, nMedian

    for iii in range(nData):
        dyd[iii] = dyv[iii]

    for iii in range(nMedian,int(nData-(nMedian+1)), int(nMedian/2)):
        if dyv[iii] > 1.*dyrms:
#            print 'range: %d to %d (%d) ' % (iii-nMedian, iii+nMedian, nMedian)
            xs = xv[(iii-nMedian):(iii+nMedian)]
            ys = dyv[(iii-nMedian):(iii+nMedian)]
            n = len(xs)
            x0 = xv[iii]
            y0 = dyv[iii]
#        sigma = dabs(sum(ys*(xs-x0)**2)/n)        #note this correction
            print 'Array sizes: ', len(xs), y0, x0, (xv[iii]-xv[iii-2])
            popt,pcov = curve_fit(gaus,xs,ys,p0=[y0, x0, (xv[iii]-xv[iii-2])])
            print iii, xv[iii], dyv[iii]
            print iii, popt, '+/-', pcov[0]
# now subtract the fit
            for jjj in range( -nMedian,nMedian):
                # compute the correction to the fit
                ady = gaus( xv[jjj+iii], popt[0], popt[1], popt[2])
                if ady > 1:
                    print iii+jjj, xv[jjj+iii], dyv[jjj+iii], ady
                dyd[jjj+iii] = dyv[jjj+iii] - ady

# transfer data with values flagged        
        for iii in range(nData):
            hv[iii] = dyd[iii]
# step through the median values
    plt.plot(xv, dyv, colors[nloop], linestyle=linestyles[nloop],label="Median-Fit")
    nloop = nloop+1
    nMedian = nMedian - 2

plt.plot(xv, medianValues, 'r', linestyle=linestyles[0],label="Median")
plt.legend(loc='upper right')
plt.show()
# exit while median loop


#print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
#plt.plot(xv, hv, colors[0], linestyle=linestyles[0],label="hot")
#plt.plot(xv, cv, colors[1], linestyle=linestyles[0],label="cold")
#plt.plot(xv, gain, colors[1], linestyle=linestyles[0],label="gain")

trx = np.zeros(nData)
for iii in range(nData):
    trx[iii] = (cv[iii]/gain[iii]) - tcold

Tsys = statistics.median(trx)
print "Median Receiver + Antenna Temp: ", Tsys


