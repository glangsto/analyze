#Python Script to plot calibrated/baseline-fit  NSF record data.
#plot the raw data from the observation
#HISTORY
#19SEP11 GIL do not use statistics, use numpy equivalent
#18DEC11 GIL minor code cleanup 
#18MAR05 GIL implement folding option
#18FEB06 GIL keep track of first az,el
#17AUG18 GIL plot the last spectrum in the series
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
import radioastronomy
import copy
import interpolate

# default values
avetimesec = 3600.
# put your list of known RFI features here.  Must have at least two.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]

nargs = len(sys.argv)
#first argument is the averaging time in seconds
timearg = 1
namearg = 2

# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False
flagCenter = True
# put list of RFI features here, for interpolation later
flagRfi = False
flagRfi = True
# to address an old problem, optionally allow folding spectra
doFold = False
# optionally turn on debug plotting

doDebug = False
iarg = 1
# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-F':
        print 'Folding specectra'
        doFold = True
    elif sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-D':
        print 'Adding Debug Printing'
        doDebug = True
    else:
        break
    iarg = iarg + 1
    timearg = timearg+1
    namearg = namearg+1
# end of while not reading file names

#first argument is the averaging time in seconds
avetimesec = float(sys.argv[timearg])
print "Average time: ", avetimesec, " (seconds)"
newObs = False

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-c','-m','-y','-c','-m','-y','-c','-m','-y','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
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
maxvel = 200.
minvel = -180.
# currently used
maxvel = 180.
minvel = -550.

c = 299792.  # (v km/sec)
nuh1 = 1420.4056 # neutral hydrogen frequency (MHz)
thot = 285.0  # define hot and cold
#thot = 272.0  # 30 Farenheit = 272 K
tcold = 10.0
tmin = 20.0 
tmax = 999.0 # define reasoanable value limits

# prepare to remove a linear baseline

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
firstaz = -1
otheraz = -1
dt = datetime.timedelta(seconds=0.)

# rest of arguments are file names
names = sys.argv[namearg:]
names = sorted(names)
nFiles = len(names)
# create the spectrum class/structure to receive spectra
rs = radioastronomy.Spectrum()

for filename in names:

    parts = filename.split('/')
    nparts = len(parts)
    if nparts == 1:
        aname = parts[0]
    else:
        aname = parts[nparts-1]
    parts = aname.split('.')
    nparts = len(parts)
    if nparts < 2:
        print 'File is not an astronomy file: ',filename
        continue
    else:
        extension = parts[nparts-1]
    extension = extension.upper()
    if (extension != 'HOT') and (extension != 'AST') and (extension != 'CLD'):
        print 'Extension not recognized : ', parts[nparts-1]
        continue

    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

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
        if minel > rs.telel:
            minel = rs.telel
        if maxel < rs.telel:
            maxel = rs.telel
        if firstaz < 0:
            firstaz = rs.telaz
            otheraz = firstaz
        if firstaz != rs.telaz:
            otheraz = rs.telaz
            
if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA / float(nhot)
    print "Found %3d Hot load observations" % nhot
else:
    print "No Hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.E-6
yv = hot.ydataA
nData = len( xv)        # get data length and indicies for middle of spectra 
n6 = int(nData/6)
n56 = 5*n6

xa = 200
#xb = 1024-xa
xb = nData-200


#previously just copy
#hv = yv
if flagRfi:
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

vel = np.zeros(nData)
# create index array
iv = np.zeros(nData)
for jjj in range (0, nData):
    vel[jjj] = c * (nuh1 - xv[jjj])/nuh1
    iv[jjj] = jjj

# velocity decreases with increasing channel #
for jjj in range (1, (nData-1)):
    if (minvel < vel[jjj]) and (minvel >= vel[jjj+1]):
        xa = jjj
    if (maxvel < vel[jjj]) and (maxvel >= vel[jjj+1]):
        xb = jjj

if doDebug:
    print 'Min Vel at channel: ',xa, minvel
    print 'Max Vel at channel: ',xb, maxvel
                                   
# will need smoothed hot load values in remaining calc
#for iii in range(1,(nData-2)):
#    hv[iii] = (yv[iii-2]+yv[iii-1]+yv[iii]+yv[iii+1]+yv[iii+2])/5.
#hv = yv

if doDebug:
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

# but if no low latitude data are available use minimum
# now average coldest data for calibration
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

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
    print "Found %3d High Galactic Latitude spectra" % (nhigh)
    yv = high.ydataA

    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

# finally compute gain on a channel by channel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)
#    gain[iii] = hv[iii]/thot

#now want gain using only hot counts, but must add in Tsys
tSys = cv/gain
tSysMiddle = np.median(tSys[n6:n56])

# for remainder of calculations only use hot counts for calibration
for iii in range(nData):
    gain[iii] = hv[iii]/(thot + tSysMiddle - tcold)

def compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold):
    nData = len(xv)
    vel = np.zeros(nData)
    tsys = np.zeros(nData)    # initialize arrays
    tsky  = np.zeros(nData)    # initialize arrays

    for jjj in range (0, nData):
        vel[jjj] = c * (nuh1 - xv[jjj])/nuh1
        tsys[jjj] = yv[jjj]/gain[jjj]
        tsky[jjj] = tsys[jjj] - Tsys

    if flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
        icenter = int(nData/2)
        tsky[icenter] = (tsky[icenter-2] + tsky[icenter+2])*.5
        tsky[icenter-1] = (3.*tsky[icenter-2] + tsky[icenter+2])*.25
        tsky[icenter+1] = (tsky[icenter-2] + 3.*tsky[icenter+2])*.25

    iref = int(nData/2)
    vref = vel[iref]
    dv   = (vel[iref+2]-vel[iref-2])/4.
    imin = int(((minvel - vref)/dv) + iref)
    imax = int(((maxvel - vref)/dv) + iref) + 1
    if imax < imin: # swap indices if frequency opposit velocity
        temp = imin
        imin = imax
        imax = temp

    if imin < 0:
        print 'Imin Error computing baseline: ', imin
        imin = 0
    if imin >= nData:
        print 'Imin Error computing baseline: ', imin
        imin = nData-1

    if imax < 0:
        print 'Imax Error computing baseline: ', imax
        imax = 0
    if imax >= nData:
        print 'Imax Error computing baseline: ', imax
        imax = nData-1

    ya = np.median(tsky[(xa-10):(xa+10)])
    yb = np.median(tsky[(xb-10):(xb+10)])
    slope = (yb-ya)/(xb-xa)
#                baseline = tsky
#                print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
    for iii in range( nData):
        tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))

    return tsky, vel

fig, ax1 = plt.subplots(figsize=(10, 6))
#plt.hold(True)
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
ymed = np.median(yv[n6:n56])
ystd = np.std(yv[n6:n56])
if ystd <= 0.:
    ystd = 0.001
count = hot.count
ncold = 0

# set indicies for normalizing intensities
bData = 790
eData = 900
bData = int(3*int(nData/10))
eData = int(6.5*int(nData/10))

# compute median raw values
hotmedian = np.median(hv[bData:eData])
coldmedian = np.median(cv[bData:eData])
hvnorm = (1./hotmedian) * hv
cvnorm = (1./coldmedian) * cv

#
deltanorm = cvnorm - hvnorm
deltanorm = deltanorm + 1.0
zeronorm = 0.0*deltanorm + 1.0

# compute the reciever temperature
trx = np.zeros(nData)
zeros = np.zeros(nData)
for iii in range(nData):
    trx[iii] = (cv[iii]/gain[iii]) - tcold

Tsys = np.median(trx[bData:eData])
print "Median Receiver + Antenna Temp: %7.2f (K)" % ( Tsys)

#plt.plot(xv, trx, colors[1], linestyle=linestyles[0],label="Tsys")
#plt.legend(loc='upper left')
#plt.show()

avetime = datetime.timedelta(seconds=avetimesec)

#plt.plot(vel, hv, colors[0], linestyle=linestyles[3],label="hot")

nRead = 0        

# now read through all data and average cold sky obs
for filename in names:

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    extension = parts[1]
    nRead = nRead + 1
    if doDebug:
        print "%3d: %s " % (nRead, filename)

# exclude hot load data for averaging
    if extension == 'hot':
        continue

    rs = radioastronomy.Spectrum()
#  print filename
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

    # recreate time/date string from UTC
    autc = str(rs.utc)
    parts = autc.split(' ')
    date = parts[0]
    nd = len(date)
    date = date[2:nd]
    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format
    parts = time.split('.')  # trim off seconds part of time
    time = parts[0]

    if firstdate == "":
        firstdate = date


# if a sky observation
    if rs.telel > 0.:

# if first time reading data, set obs parameters
        if lastfreq == 0.:
            lastfreq = rs.centerFreqHz 
            lastbw = rs.bandwidthHz
            lastgain = rs.gains[0]
            lastaz = rs.telaz
            lastel = rs.telel
            cold = copy.deepcopy( rs)
            crosszero = False
            firstlon = rs.gallon
            ncold = 0

        if ncold > 0:
            # time difference is between mid-points of integrations
            dt = rs.utc - cold.utc 
            # add the time since midpoint of latests
            dt = dt + datetime.timedelta(seconds=rs.durationSec/2.)
            # plus time before start of the first
            dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)

            lastdate = date

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

        if nRead == nFiles:   # if this is the last file, must force output
            newObs = True

        # if time to average (or end of all files)
        if (dt > avetime) or (filename == sys.argv[nargs-1]) or newObs:
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
            if flagRfi:
                yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

            xmin = min(xv)
            xmax = max(xv)
            xallmin = min(xmin, xallmin)
            xallmax = max(xmax, xallmax)
            count = cold.count
            note = cold.noteA
                    #print('%s' % note)
            ncolor = min(nmax-1, nplot) 

            label = 'L,L=%5.1f,%5.1f' % (gallon, gallat)
            
            tsky, vel = compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold)
            ymed = np.median(tsky[n6:n56])
            ystd = np.std(tsky[n6:n56])
            if ystd <= 0.:
                ystd = 0.001
            
            ymin = min(tsky[(1*nData/4):(3*nData/4)])
            ymax = max(tsky[(1*nData/4):(3*nData/4)])
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            avedatetime = cold.utc.isoformat()
            dates = avedatetime.split('T')
            atime = dates[1]
            timeparts = atime.split('.')
            labeltime = timeparts[0]
            label = '%s, A,E: %5s,%5s, L,L: %5.1f,%5.1f' % (labeltime, az, el, gallon, gallat)
            if minel == maxel:
                label = '%s L,L=%5.1f,%5.1f' % (labeltime, gallon, gallat)
            else:
                label = '%s L,L=%5.1f,%4.1f A,E=%4.0f,%4.0f' % (labeltime, gallon, gallat, az, el)
            if ncold > 0:
                print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, (ymax-ymed)/ystd, ncold, label)
                if gallat < 7.5 and gallat > -7.5:
                    plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
                elif gallat < 15. and gallat > -15.:
                    plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
                else:
                    plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
                nplot = nplot + 1
            ncold = 0

            # end if a cold file
            if nRead > nFiles:
                break

# if this was a new obs; restart the sums
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
            timesum = rs.durationSec
        else: # else not enough time yet, average cold data
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
            timesum = timesum + rs.durationSec
#                plt.plot(iss, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            # end if not a enough time

        # end if a cold file
        if nRead > nFiles:
            break

    #end for all files to sum
# end of all files to read

if doDebug:
    print 'Number of remaining observations not plotted: ', ncold

# if data remaing from summation
if ncold > 0:
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
    if flagRfi:
        yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
    
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin, xallmin)
    xallmax = max(xmax, xallmax)
    count = cold.count
    note = cold.noteA
    ncolor = min(nmax-1, nplot) 

    tsky, vel = compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold)
    ymed = np.median(tsky[n6:n56])
    ystd = np.std(tsky[n6:n56])
    if ystd <= 0.:
        ystd = 0.001
            
    ymin = min(tsky[n6:n56])
    ymax = max(tsky[n6:n56])
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    avedatetime = cold.utc.isoformat()
    dates = avedatetime.split('T')
    atime = dates[1]
    timeparts = atime.split('.')
    labeltime = timeparts[0]
    if minel == maxel:
        label = '%s L,L=%5.1f,%5.1f' % (labeltime, gallon, gallat)
    else:
        label = '%s L,L=%5.1f,%5.1f A,E=%4.0f,%4.0f' % (labeltime, gallon, gallat, az, el)
    print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, (ymax-ymed)/ystd, ncold, label)
    if gallat < 7.5 and gallat > -7.5:
        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
    elif gallat < 15. and gallat > -15.:
        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
    else:
        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
    nplot = nplot + 1
        
# if observations cover multiple days
if firstdate != lastdate:
    date = firstdate + " - " + lastdate
else:
    date = firstdate

mytitle = "%s    " % (date)
if (firstaz == otheraz):
    mytitle = mytitle + "Az = %6.1f, " % (firstaz)
else:
    mytitle = mytitle + "Az = %6.1f to %6.1f, " % (firstaz, otheraz)

if minel == maxel:
    mytitle = mytitle + " El=%6.1f" % (minel)
else:
    mytitle = mytitle + " El=%6.1f to %6.1f" % (minel, maxel)

fig.canvas.set_window_title(mytitle)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
#plt.xlim(xallmin,xallmax)  # normal plot range
#plt.xlim(-150., 300.)
# select the relevant range of velocities (km/sec) for plotting
#plt.xlim(-400., 250.)      # set velocity range for particular investigation
plt.xlim(-600., 300.)
# keep the plot from becoming too narrow
if yallmax < 8:
    yallmax = 8
# set the y scale to go above and below all data
plt.ylim((yallmin*.97)-1., 1.15*yallmax)
plt.title(mytitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.show()
