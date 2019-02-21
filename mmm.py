#Python Script to comput and plot the min, max, median and average raw NSF record data.
#HISTORY
#17SEP22 GIL init version based on m.py
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import statistics
import radioastronomy
import copy
from scipy.signal import savgol_filter
import interpolate

dy = -1.
linelist = [1420.0, 1419.0, 1418.0]  # RFI lines in MHz
linewidth = [9, 5, 7]

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
firstfile = 1

# first read through all data and find hot load
names = sys.argv[firstfile:]
names = sorted(names)
nfiles = len(names)
nData = 0
nRead = 0
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    
    # if spectra array is not yet initialized
    if nData == 0:
        nData = len(rs.xdata)
        spectra = np.zeros((nfiles,nData))
        rsmin = copy.deepcopy( rs)
        rsmax = copy.deepcopy( rs)
        rsmed = copy.deepcopy( rs)
        rsave = copy.deepcopy( rs)

# convert to MHz
    xv = rs.xdata*1.E-6
    yv = rs.ydataA * scalefactor
    if doFlag:
        hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
        spectra[nRead, : ] = hv
    else:
        spectra[nRead, : ] = yv
    nRead = nRead + 1

print 'Read ',nRead,' Spectra with ',nData,' Spectral channels '

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
    
fig, ax1 = plt.subplots(figsize=(10, 6))
#plt.hold(True)

plt.plot(xv, rsmax.ydataA, colors[1], linestyle=linestyles[1],label='Maximum', lw=4)
plt.plot(xv, rsmed.ydataA, colors[2], linestyle=linestyles[2],label='Median', lw=4)
plt.plot(xv, rsmin.ydataA, colors[0], linestyle=linestyles[0],label='Minimum', lw=4)
plt.plot(xv, rsave.ydataA, colors[3], linestyle=linestyles[3],label='Average', lw=4)

print 'Maximum Spectrum File:'
print names[imax]

rs.read_spec_ast(names[imax])
gallon = rs.gallon
gallat = rs.gallat
label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( names[imax],rs.telaz,rs.telel,gallon,gallat)
plt.plot(xv, rs.ydataA*scalefactor, colors[5], linestyle=linestyles[5],label=label, lw=4)
rs.read_spec_ast(names[imin])
gallon = rs.gallon
gallat = rs.gallat
label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( names[imin],rs.telaz,rs.telel,gallon,gallat)
plt.plot(xv, rs.ydataA*scalefactor, colors[4], linestyle=linestyles[4],label=label, lw=4)
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
