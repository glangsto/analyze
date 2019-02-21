#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#17AUG17 GIL Note elevation range
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
import interpolate

avetimesec = 3600.

avetimesec = 900.
#avetimesec = 60.
avetimesec = 120.
dy = -1.


linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

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
lastel = 0.
lastaz = 0.
firstdate = ""
lastdate = ""
minel = 200.
maxel = -200.

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
    else:
        if minel > rs.telel:
            minel = rs.telel
        if maxel < rs.telel:
            maxel = rs.telel

if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA / float(nhot)
    print "Found %3d hot load obs" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.E-6
yv = hot.ydataA
#yv = hot.foldfrequency()
nData = len(xv)

vel = np.zeros(nData)
for jjj in range (0, nData):
    vel[jjj] = c * (xv[jjj] - nuh1)/nuh1

hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
yv = copy.deepcopy(hv)

#interpolate to smooth reference 
#for iii in range(1,(nData-2)):
#    hv[iii] = (yv[iii-2]+yv[iii-1]+yv[iii]+yv[iii+1]+yv[iii+2])/5.

fig, ax1 = plt.subplots(figsize=(10, 6))
#plt.hold(True)
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
ymed = statistics.median(yv)
count = hot.count
ncold = 0

# condition data for avoid divide by zero later (2 depends on scalefactor)
for iii in range(nData):
    if hv[iii] < 2.:
        hv[iii] = 2.

avetime = datetime.timedelta(seconds=avetimesec)

nread = 0        
nFiles = len(names)
print nFiles

# now read through all data and average cold sky obs
for filename in names:

    rs = radioastronomy.Spectrum()
#    print filename
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    parts = aname.split('T')
    date = parts[0]
    if firstdate == "":
        firstdate = date
    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format
    nread = nread + 1

# if a sky observation
    if rs.telel > 0.:

# if first time reading data, set obs parameters
        if lastfreq == 0.:
            lastfreq = rs.centerFreqHz 
            lastbw = rs.bandwidthHz
            lastgain = rs.gains[0]
            cold = copy.deepcopy( rs)
            ncold = 0

        if ncold > 0:
            # time difference is between mid-points of integrations
            dt = rs.utc - cold.utc 
            # add the time since midpoint of latests
            dt = dt + datetime.timedelta(seconds=rs.durationSec/2.)
            # plus time before start of the first
            dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)

            newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
            newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl
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
                if newAzEl:
                    print "LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees"
                    lastaz = rs.telaz
                    lastel = rs.telel

            print 'Nread, nFiles: ', nread, nFiles
            if nread == nFiles:   # if this is the last file, must force output
                newObs = True
                    
           # if time to average (or end of all files)
            if (dt > avetime) or (filename == sys.argv[nargs-1]) or newObs:
                cold.ydataA = cold.ydataA/float(timesum)

                gallon = cold.gallon/float(timesum)
                gallat = cold.gallat/float(timesum)
                az = cold.telaz
                el = cold.telel

                # convert to MHz from Hz
                xv = cold.xdata * 1.E-6
                yv = cold.ydataA * scalefactor

                yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
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

                ymin = min(tsky[(nData/8):(7*nData/8)])
                ymax = max(tsky[(nData/8):(7*nData/8)])
                yallmin = min(ymin,yallmin)
                yallmax = max(ymax,yallmax)
                ymed = statistics.median(tsky)
                if minel == maxel: 
                    label = '%s Lon,Lat=%5.1f,%5.1f' % (time,gallon,gallat)
                else:
                    label = '%s Lon,Lat=%5.1f,%5.1f El=%4.1f' % (time, gallon, gallat, el)
                print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)

                # plot thicker lines when near the galactic plane
                if gallat < 7.5 and gallat > -7.5:
                    plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
                elif gallat < 15. and gallat > -15.:
                    plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
                else:
                    plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)

                nplot = nplot + 1
                ncold = 0
            # end if a new observation
# if this was a new obs; restart the sums
        if ncold == 0:
            cold = rs  # initial spectrum is one just read
            ncold = 1
#            print 'Xmin: ', min(cold.xdata)/1e6, 'Xmax: ', max(cold.xdata),' MHz'
            # sums are weighted by durations
            firstlon = rs.gallon
            cold.ydataA = cold.ydataA * cold.durationSec
            cold.gallat = cold.gallat * cold.durationSec
            cold.gallon = cold.gallon * cold.durationSec
            # keep track of observing time for weighted sum
            timesum = cold.durationSec
        else: # else ont enough time yet, average cold data
            cold.count = cold.count + rs.count
            ncold = ncold + 1
            cold.ydataA = cold.ydataA + (rs.ydataA * cold.durationSec)
            # fix wrap of longitudes
            if abs(rs.gallon - firstlon) > 180:
                crosZero = True
                if rs.gallon > firstlon:
                    rs.gallon = rs.gallon - 360.
                else:
                    rs.gallon = rs.gallon + 360.
            cold.gallon = cold.gallon + (rs.gallon * cold.durationSec)
            cold.gallat = cold.gallat + (rs.gallat * cold.durationSec)
            # keep track of observing time for weighted sum
            timesum = timesum + cold.durationSec
            # end if not a enough time
        # end if a cold file
    #end for all files to sum

#plt.xlim(xallmin,xallmax)
lastdate = date
# if observations span several days
if lastdate != firstdate:
    date = firstdate + "-" + lastdate

# if change in elevation 
if minel == maxel:
    mytitle = "%s    Az=%6.1f, El=%6.1f" % (date, az, minel)
else:
    mytitle = "%s    Az=%6.1f, El=%6.1f to %6.1f" % (date, az, minel, maxel)

fig.canvas.set_window_title(mytitle)

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
#plt.xlim(-300., 500.)
plt.xlim(xallmin, xallmax)
plt.ylim((yallmin*.97)-1., 1.25*yallmax)
#plt.ylim(yallmin*.97, 400.)
plt.title(mytitle, fontsize=16)
plt.xlabel('Frequency (MHz)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
