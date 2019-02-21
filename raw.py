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
#import radioastronomy
import tsys
import copy

avetimesec = 120.

nargs = len(sys.argv)

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

tcal = tsys.Tsys()

# first read through all data and find hot load
for iii in range(2, nargs):

    filename = sys.argv[iii]

    if filename == "Off.ast":
        tcal.off.read_spec_ast(filename)
        tcal.off.azel2radec()    # compute ra,dec from az,el
        print 'Read an Off file'
        continue
    else:
        tcal.hot.read_spec_ast(filename)
        tcal.hot.azel2radec()    # compute ra,dec from az,el

    
    if tcal.hot.telel < 0:
       if nhot == 0:
            hot = copy.deepcopy(tcal.hot)
#            hot.ydataA = tcal.hot.ydataA
            nhot = 1
       else:
            hot.ydataA = hot.ydataA + tcal.hot.ydataA
            hot.count = hot.count + tcal.hot.count
            nhot = nhot + 1

if nhot > 0:
    tcal.hot.ydataA = scalefactor * hot.ydataA / float(nhot)
    print "Found %3d hot load obs" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

# convert to MHz
xv = tcal.hot.xdata * 1.E-6
yv = tcal.hot.ydataA
nData = len(xv)
hv = yv

#for iii in range(2,(nData-2)):
#    hv[iii] = (yv[iii-2]+yv[iii-1]+yv[iii]+yv[iii+1]+yv[iii+2])/5.

#smooth hot slightly
for iii in range(1,(nData-1)):
    hv[iii] = (yv[iii-1]+yv[iii]+yv[iii+1])/3.

fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
az = tcal.hot.telaz
el = tcal.hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
ymed = statistics.median(yv)
count = hot.count
ncold = 0
#print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
#plt.plot(xv, yv, colors[0], linestyle=linestyles[0],label=label)

vel = np.zeros(nData)
for jjj in range (0, nData):
    vel[jjj] = c * (xv[jjj] - nuh1)/nuh1

# condition data for avoid divide by zero later (2 depends on scalefactor)
for iii in range(nData):
    if hv[iii] < 2.:
        hv[iii] = 2.

plt.plot(vel, hv, colors[0], linestyle=linestyles[3], label="hot")

avetime = datetime.timedelta(seconds=avetimesec)

nread = 0        
# now read through all data and average cold sky obs
for iii in range(2, nargs):

    filename = str(sys.argv[iii])
    if (filename == "Off.ast"):
        continue
    
#    print filename
    tcal.cold.read_spec_ast(filename)
    tcal.cold.azel2radec()    # compute ra,dec from az,el

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    parts = aname.split('T')
    date = parts[0]
    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format

# if a sky observation
    if tcal.cold.telel > 0.:

# if first time reading data, set obs parameters
        if lastfreq == 0.:
            lastfreq = tcal.cold.centerFreqHz 
            lastbw = tcal.cold.bandwidthHz
            lastgain = tcal.cold.gains[0]
            timesum = tcal.cold.durationSec
            cold = copy.deepcopy( tcal.cold)
            ncold = 0

        if ncold > 1:
            # time difference is between mid-points of integrations
            dt = tcal.cold.utc - cold.utc 
            # add the time since midpoint of latests
            dt = dt + datetime.timedelta(seconds=tcal.cold.durationSec/2.)
            # plus time before start of the first
            dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)

            newObs = (lastfreq != tcal.cold.centerFreqHz) or (lastbw != tcal.cold.bandwidthHz) or (lastgain != tcal.cold.gains[0])
            if newObs:
                print "Change in observing parameters: "
                if lastfreq != tcal.cold.centerFreqHz:
                    print "LastFreq: ", lastfreq/1e6, "New: ", tcal.cold.centerFreqHz/1e6, " MHz"
                    lastfreq = tcal.cold.centerFreqHz
                if lastbw != tcal.cold.bandwidthHz:
                    print "LastBandwidth: ", lastbw/1e6, "New: ", tcal.cold.bandwidthHz/1e6, " MHz"
                    lastbw = tcal.cold.bandwidthHz
                if lastgain != tcal.cold.gains[0]:
                    print "LastGain: ", lastgain, "New: ", tcal.cold.gains[0], " dB"
                    lastgain = tcal.cold.gains[0]
           # if time to average (or end of all files)
            if (dt > avetime) or (iii == (nargs-1)) or newObs:
                cold.ydataA = cold.ydataA/float(timesum)
                cold.gallon = cold.gallon/float(timesum)
                cold.gallat = cold.gallat/float(timesum)
                az = cold.telaz
                el = cold.telel

                tcal.cold = copy.deepcopy(cold)

                # convert to MHz
                xv = cold.xdata * 1.E-6
                yv = cold.ydataA * scalefactor

                xmin = min(xv)
                xmax = max(xv)
                xallmin = min(xmin, xallmin)
                xallmax = max(xmax, xallmax)
                count = cold.count
                note = cold.noteA
                    #print('%s' % note)
                ncolor = min(nmax-1, nplot) 

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

                label = 'Lon,Lat=%5.1f,%5.1f' % (tcal.cold.gallon, tcal.cold.gallat)

                ymin = min(tsky[(nData/8):(7*nData/8)])
                ymax = max(tsky[(nData/8):(7*nData/8)])
                yallmin = min(ymin,yallmin)
                yallmax = max(ymax,yallmax)
                ymed = statistics.median(tsky)
                label = '%s, Az,El: %5s,%5s, Lon,Lat: %5.1f,%5.1f' % (time, az, el, tcal.cold.gallon, tcal.cold.gallat)
                label = '%s Lon,Lat=%5.1f,%5.1f' % (time, tcal.cold.gallon, tcal.cold.gallat)
                print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)
#                plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
#                print 'Y min: ', min(yv), 'Y max :', max(yv)
                plt.plot(vel, yv, colors[ncolor], linestyle=linestyles[ncolor],label=label)

                nplot = nplot + 1
                ncold = 0
            # end if a new observation
# if this was a new obs; restart the sums
        if ncold == 0:
            cold = copy.deepcopy(tcal.cold)  # initial spectrum is one just read
            cold.xdata = tcal.cold.xdata
            print 'Xmin: ', min(cold.xdata)/1e6, 'Xmax: ', max(cold.xdata),' MHz'
            print 'Xmin: ', min(tcal.cold.xdata)/1e6, 'Xmax: ', max(tcal.cold.xdata),' MHz'
            ncold = 1
            # sums are weighted by durations
            cold.ydataA = tcal.cold.ydataA * tcal.cold.durationSec
            cold.gallat = tcal.cold.gallat * tcal.cold.durationSec
            cold.gallon = tcal.cold.gallon * tcal.cold.durationSec
            # keep track of observing time for weighted sum
            timesum = tcal.cold.durationSec
        else: # else ont enough time yet, average cold data
            cold.count = cold.count + tcal.cold.count
            ncold = ncold + 1
            cold.xdata = tcal.cold.xdata
            cold.ydataA = cold.ydataA + (tcal.cold.ydataA * tcal.cold.durationSec)
            cold.gallon = cold.gallon + (tcal.cold.gallon * tcal.cold.durationSec)
            cold.gallat = cold.gallat + (tcal.cold.gallat * tcal.cold.durationSec)
            # keep track of observing time for weighted sum
            timesum = timesum + tcal.cold.durationSec
#                plt.plot(iss, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            # end if not a enough time
        # end if a cold file
    #end for all files to sum

#plt.xlim(xallmin,xallmax)
fig.canvas.set_window_title(date)
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
#plt.xlim(-50., 800.)
plt.ylim(yallmin*.97, 2*yallmax)
#plt.ylim(0, max(hv)*1.2)
plt.title(date, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.show()
