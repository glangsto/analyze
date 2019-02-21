#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#16Oct07 GIL check for changes in frequency and/or gain
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import statistics
import radioastronomy
import copy
#import list

lonavebegin = 350.
lonaveend   = 5.

avetimesec = 900.
#avetimesec = 60.
avetimesec = 120.
dy = -1.

nargs = len(sys.argv)
names = sys.argv[6:]
names = sorted(names)

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1e8
xallmax = -9.e9
xallmin = 9.e9
yallmax = -9.e9
yallmin = 9.e9

c = 299792.  # (v km/sec)
nuh1 = 1420.4056 # neutral hydrogen frequency (MHz)
thot = 285.0  # define hot and cold
#thot = 272.0  # 30 Farenheit = 272 K
tcold = 10.0
tmin = 20.0 
tmax = 999.0 # define reasoanable value limits

# prepare to remove a linear baseline
xa = 100
xb = 1024-100

#for symbol, value in locals().items():
#    print symbol, value

nplot = 0
nhot = 0
ncold = 0
lastfreq = 0.
lastbw = 0.
lastgain = 0.

#first argument is the averaging time in seconds

avetimesec = float(sys.argv[1])
print "Average time: ", avetimesec, " (seconds)"
LABEL = 0
lonavebegin = float(sys.argv[2])
lonaveend = float(sys.argv[3])
latavebegin = float(sys.argv[4])
lataveend = float(sys.argv[5])
print "Average longitude; %7.1f to %7.1f (degrees)" % (lonavebegin, lonaveend)
print "Average latitude;  %7.1f to %7.1f (degrees)" % (latavebegin, lataveend)
label = ["",""]
label[LABEL] = "plot label"

def plotcold( time, timesum, xa, xb, xallmin, xallmax, cold):
    """
    Compute the plot based on the averaged data
    """
    cold.ydataA = cold.ydataA/float(timesum)

    gallon = cold.gallon/float(timesum)
    if gallon > 360.:
        gallon = gallon - 360.
    elif gallon < 0.:
        gallon = gallon + 360.
        
    gallat = cold.gallat/float(timesum)
    az = cold.telaz
    el = cold.telel
    
    xv = cold.xdata * 1.E-6
    yv = cold.ydataA * scalefactor

    nData = len(xv)
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin, xallmin)
    xallmax = max(xmax, xallmax)
    count = cold.count
    note = cold.noteA

    tsys = np.zeros(nData)    # initialize arrays
    iss = np.zeros(nData)
    Z = np.zeros(nData)
    oneMZ = np.zeros(nData)
    for jjj in range(0, nData):
        iss[jjj] = jjj
        if yv[jjj] < 1.:
            yv[jjj] = 1.
            
        Z[jjj] = yv[jjj]/hv[jjj]
        oneMZ[jjj] = 1. - Z[jjj]
        if oneMZ[jjj] < .001:
            oneMZ[jjj] = 0.001
        tsys[jjj] = ((Z[jjj]*thot) - tcold)/oneMZ[jjj]

    vel = np.zeros(nData)
    for jjj in range (0, nData):
        vel[jjj] = c * (xv[jjj] - nuh1)/nuh1

    tsky  = np.zeros(nData)    # initialize arrays
    S     = np.zeros(nData)    # initialize arrays

    for jjj in range (0, nData):
        vel[jjj] = c * (xv[jjj] - nuh1)/nuh1
        S[jjj] = hv[jjj]/(tsys[jjj] + thot)
        tsky[jjj] = yv[jjj]/S[jjj]

        # remove spike in center of the plot
    icenter = 512
    tsky[icenter] = (tsky[icenter-2] + tsky[icenter+2])*.5
    tsky[icenter-1] = (3.*tsky[icenter-2] + tsky[icenter+2])*.25
    tsky[icenter+1] = (tsky[icenter-2] + 3.*tsky[icenter+2])*.25

    label = 'Lon,Lat=%5.1f,%5.1f' % (gallon, gallat)

    ymed = statistics.median(tsky)
    ya = statistics.median(tsky[(xa-10):(xa+10)])
    yb = statistics.median(tsky[(xb-10):(xb+10)])
    slope = (yb-ya)/(xb-xa)
#        print 'slope: ', slope
    for iii in range( 0, nData):
        tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))

    label = '%s, Az,El: %5s,%5s, Lon,Lat: %5.1f,%5.1f' % (time, az, el, gallon, gallat)
    label = '%s Lon,Lat=%5.1f,%5.1f' % (time,gallon,gallat)
    print '%d: Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ncold, ymax, ymed, ymax/ymed, ncold, label)
    cold.xdata = vel
    cold.ydataA = tsky
    print 'Label: ',label
    return label
        
# first read through all data and find hot load
for iii in range(6, nargs):

    filename = sys.argv[iii]
    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    extension = parts[1]
# now only processing hot loads
    if extension != 'hot':
#        print parts[1], extension
        continue
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

if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA / float(nhot)
    print "Found %3d hot load obs" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.E-6
yv = hot.ydataA
nData = len(xv)
hv = yv

vel = np.zeros(nData)
for jjj in range (0, nData):
    vel[jjj] = c * (xv[jjj] - nuh1)/nuh1

# will need smoothed hot load values in remaining calc
#for iii in range(1,(nData-2)):
#    hv[iii] = (yv[iii-2]+yv[iii-1]+yv[iii]+yv[iii+1]+yv[iii+2])/5.
hv = yv

fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
ymed = statistics.median(yv)
count = hot.count
ncold = 0
#print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
#plt.plot(xv, yv, colors[0], linestyle=linestyles[0],label=label)

# condition data for avoid divide by zero later (2 depends on scalefactor)
for iii in range(nData):
    if hv[iii] < 2.:
        hv[iii] = 2.

avetime = datetime.timedelta(seconds=avetimesec)
time = ""
date = ""

#plt.plot(vel, hv, colors[0], linestyle=linestyles[3],label="hot")

nread = 0        

# now read through all data and average cold sky obs
for aname in names:

    filename = str(aname)

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    parts = aname.split('T')
    extension = parts[1]
    if extension == "hot":
        continue

    rs = radioastronomy.Spectrum()
#  print filename
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el

# if a sky observation
    if rs.telel > 0.:

        el = rs.telel
        az = rs.telaz

# if longitude is in range for average
        # if average does not cross zero
        if lonaveend > lonavebegin:
            if (rs.gallon > lonaveend) or (rs.gallon < lonavebegin):
#                print 'A Skipping ',iii, rs.gallon,rs.gallat 
                continue
        else: # else average crosses zero
            if (rs.gallon < lonaveend) and (rs.gallon > lonavebegin):
#                print 'B Skipping ',iii, rs.gallon,rs.gallat 
                continue
        # if latitide is in range
        if (rs.gallat > lataveend) or (rs.gallat < latavebegin):
#            print 'C Skipping ',iii, rs.gallon,rs.gallat 
            continue

# if first time reading data, set obs parameters
        if lastfreq == 0.:
            lastfreq = rs.centerFreqHz 
            lastbw = rs.bandwidthHz
            lastgain = rs.gains[0]
            cold = copy.deepcopy( rs)
            crosszero = False
            firstlon = rs.gallon
            ncold = 0
            date = parts[0]
            time = parts[1]
            time = time.replace('_', ':')  # put time back in normal hh:mm:ss format


        if ncold > 1:
            # time difference is between mid-points of integrations
            dt = rs.utc - cold.utc 
            # add the time since midpoint of latests
            dt = dt + datetime.timedelta(seconds=rs.durationSec/2.)
            # plus time before start of the first
            dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)

            newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0])
            if newObs:
                print "Change in observing parameters: "
                if lastfreq != rs.centerFreqHz:
                    print "LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz"
                    lastfreq = rs.centerFreqHz
                if lastbw != rs.bandwidthHz:
                    print "LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz"
                    lastbw = rs.bandwidthHz
                if lastgain != rs.gains[0]:
                    print "LastGain: ", lastgain, "New: ", rs.gains[0], " dB"
                    lastgain = rs.gains[0]
           # if time to average (or end of all files)
            if (dt > avetime) or newObs:
                label[LABEL] = plotcold( date+" "+time, timesum, xa, xb, xallmin, xallmax, cold)
                ymin = min(cold.ydataA[int(nData/8):int(7*nData/8)])
                ymax = max(cold.ydataA[int(nData/8):int(7*nData/8)])
                yallmin = min(ymin,yallmin)
                yallmax = max(ymax,yallmax)
                    #print('%s' % note)
#                print "the values are:", label[LABEL]

                ncolor = min(nmax-1, nplot) 
                plt.plot(cold.xdata, cold.ydataA, colors[ncolor], linestyle=linestyles[ncolor],label=label[LABEL])
                nplot = nplot + 1
                ncold = 0
            # end if a new observation
# if this was a new obs; restart the sums
        if ncold == 0:
            cold = copy.deepcopy(rs)  # initial spectrum is one just read
            ncold = 1
#            print 'Xmin: ', min(cold.xdata)/1e6, 'Xmax: ', max(cold.xdata),' MHz'
            # sums are weighted by durations
            crossZero = False
            firstlon = rs.gallon
            cold.ydataA = rs.ydataA * cold.durationSec
            cold.gallat = rs.gallat * cold.durationSec
            cold.gallon = rs.gallon * cold.durationSec
            # keep track of observing time for weighted sum
            timesum = cold.durationSec
            date = parts[0]
            time = parts[1]
        else: # else ont enough time yet, average cold data
            cold.count = cold.count + rs.count
            ncold = ncold + 1
            cold.ydataA = cold.ydataA + (rs.ydataA * cold.durationSec)
            # fix wrap of longitudes
            if abs(rs.gallon - firstlon) > 180:
                crossZero = True
                if rs.gallon > firstlon:
                    rs.gallon = rs.gallon - 360.
                else:
                    rs.gallon = rs.gallon + 360.
            cold.gallon = cold.gallon + (rs.gallon * cold.durationSec)
            cold.gallat = cold.gallat + (rs.gallat * cold.durationSec)
            # keep track of observing time for weighted sum
            timesum = timesum + cold.durationSec
#                plt.plot(iss, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            # end if not a enough time
        # end if a cold file
    #end for all files to sum

if ncold > 0:
    label[LABEL] = plotcold( date+" "+time, timesum, xa, xb, xallmin, xallmax, cold)
    ymin = min(cold.ydataA[int(nData/8):int(7*nData/8)])
    ymax = max(cold.ydataA[int(nData/8):int(7*nData/8)])
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    ncolor = min(nmax-1, nplot) 
    plt.plot(cold.xdata, cold.ydataA, colors[ncolor], linestyle=linestyles[ncolor],label=label[LABEL])
    

#plt.xlim(xallmin,xallmax)
mytitle = "%s    Az=%6.1f, El=%6.1f" % (date, az, el)
fig.canvas.set_window_title(mytitle)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
plt.xlim(-250., 300.)
plt.ylim(yallmin*.97, 1.15*yallmax)
#plt.ylim( 0, max(hv))
plt.title(mytitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.show()
