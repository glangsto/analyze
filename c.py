#Python Script to plot Calibrated NSF Horn Observations
#HISTORY
#18DEC10 GIL initial version based on m.py
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import statistics
import radioastronomy
import hotcold
import copy
import interpolate

# default values
avetimesec = 3600.
# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False
flagCenter = True
# put list of RFI features here, for interpolation later
flagRfi = False
flagRfi = True
# put your list of known RFI features here.  Must have at least two.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]   # integer number of channels to interpolate over

nargs = len(sys.argv)
#first argument is the averaging time in seconds
timearg = 1
namearg = 2

doDebug = False

# if folding data
if sys.argv[timearg] == '-f':
    print 'Folding specectra'
    doFold = True
    timearg = timearg+1
    namearg = namearg+1
else:
    doFold = False
# if folding data

# if  baseline subtraciton
if sys.argv[timearg] == '-b':
    print 'Baseline subtraction'
    dosub = True
    timearg = timearg+1
    namearg = namearg+1
else:
    dosub = False

avetimesec = float(sys.argv[timearg])
print "Average time: ", avetimesec, " (seconds)"
newObs = False

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-c','-m','-y','-c','-m','-y','-c','-m','-y','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1.
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

clight = 299792.  # (v km/sec)
nuh1 = 1420.40557177667 # neutral hydrogen frequency (MHz)
nureference = 1.E6*nuh1 # neutral hydrogen frequency (Hz)
thot = 285.0  # define hot and cold in Kelvins
tcold = 10.0
tmin = 20.0  # Tsys never less than 20 K
tmax = 999.0 # define reasoanable value limits

xa = 200
xb = 1024-200
velbaseline = [-500., 200.] # array of velocities in km/sec

nplot = 0
nhot = 0         # number of obs with el < 0
ncold = 0
nhigh = 0        # highest galactic latitude
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

#first argument is the averaging time in seconds
names = sys.argv[namearg:]
names = sorted(names)
nFiles = len(names)

nhot, hot = hotcold.hotaverage( names)
if doFold:
    hot.foldfrequency()

if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA 
    print "Found %3d hot load observations" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.E-6
yv = hot.ydataA
nData = len(xv)
# if flagging RFI
if flagRfi:
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
else:
    hv = yv

if flagCenter:
    hv = hotcold.flagCenter( hv)

if doDebug:
    print "hv: ",hv[300]
vel = hot.velocities(nureference)
chanbaseline = hot.vel2chan( velbaseline, nureference)
xa = chanbaseline[0]
xb = chanbaseline[1]

# velocity decreases with increasing channel #
for jjj in range (1, (nData-1)):
    if (minvel < vel[jjj]) and (minvel >= vel[jjj+1]):
        xa = jjj
    if (maxvel < vel[jjj]) and (maxvel >= vel[jjj+1]):
        xb = jjj

print 'Min Vel at channel: ',xa, minvel
print 'Max Vel at channel: ',xb, maxvel
                                   
ncold, cold = hotcold.coldaverage( names)

if ncold < 1.:
    print "No high galactic latitude data: can not calibrate"
    exit()
else:
    cold.ydataA = scalefactor * cold.ydataA
    print "Found %d High Galactic Latitude spectra" % (ncold)
    yv = cold.ydataA

if doFold:
    cold.foldfrequency()

if flagRfi:
    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
else:
    cv = yv

if flagCenter:
    cv = hotcold.flagCenter( cv)

if doDebug:
    print "cv: ",cv[300]

fig, ax1 = plt.subplots(figsize=(10, 6))
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
ymed = statistics.median(yv)
count = hot.count
ncold = 0

# set indicies for normalizing intensities
bData = int(3*int(nData/10))
eData = int(6.5*int(nData/10))

# compute median raw values
hotmedian = statistics.median(hv[bData:eData])
coldmedian = statistics.median(cv[bData:eData])
hvnorm = (1./hotmedian) * hv
cvnorm = (1./coldmedian) * cv

#
deltanorm = cvnorm - hvnorm
deltanorm = deltanorm + 1.0
zeronorm = 0.0*deltanorm + 1.0

# compute gain in Kelvins/count
Trx, gain = hotcold.compute_gain( hv, cv, thot, tcold)


if doDebug:
    print "gv: ",gain[300]
    print "Median Receiver + Antenna Temp: ", Trx

avetime = datetime.timedelta(seconds=avetimesec)
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
# exclude hot load data for averaging
    if extension == 'hot':
        continue

    rs = radioastronomy.Spectrum()
#  print filename
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

    autc = str(rs.utc)
    parts = autc.split(' ')
    date = parts[0]
    nd = len(date)
    date = date[2:nd]
    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format
    parts = time.split('.')  # trim off seconds part of time
    time = parts[0]

# if a sky observation
    if rs.telel > 0.:
        if firstdate == "":            # if very first spectrum
            firstdate = date
            minel = rs.telel
            maxel = rs.telel
            firstaz = rs.telaz
            lastaz = rs.telaz
            otheraz = rs.telaz


        if rs.telaz != otheraz:
            otheraz = rs.telaz

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
            if flagCenter:
                yv = hotcold.flagCenter( yv)

            xmin = min(xv)
            xmax = max(xv)
            xallmin = min(xmin, xallmin)
            xallmax = max(xmax, xallmax)
            count = cold.count
            note = cold.noteA
            ncolor = min(nmax-1, nplot) 

            label = 'L,L=%5.1f,%5.1f' % (gallon, gallat)
            
            tsky = hotcold.compute_tsky( xv, yv, gain, nureference)
            ymed = statistics.median(tsky[(nData/4):(3*nData/4)])
            
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
                print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)
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

print 'Number of remaining observations not plotted: ', ncold

# if data remaing from summation
if ncold > 1:
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

    tsky, vel = hotcold.compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold, nureference)
    ymed = statistics.median(tsky)

    ymin = min(tsky[(2*nData/4):(3*nData/4)])
    ymax = max(tsky[(2*nData/4):(3*nData/4)])
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
#    print ' Max: %9.1f  Median: %9.1f Tsys: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)
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
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
# select the relevant range of velocities (km/sec) for plotting
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
