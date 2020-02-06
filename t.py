#Python Script to plot calibrated  NSF spectral integration data.
#plot the raw data from the observation
#HISTORY
#19DEC30 GIL add title option
#19SEP23 GIL use function for averaging 
#19SEP21 GIL fix finding extra spectrum
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
# specify lowest  elevation for cold load averge
lowel = 60.
# optionally turn on debug plotting
doDebug = False
myTitle = ""      # Default no input title

iarg = 1
# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-F':
        print( 'Folding specectra')
        doFold = True
    elif sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-T':   # if plot title provided
        iarg = iarg+1
        myTitle = sys.argv[iarg]
        print( 'Plot Title : ', myTitle)
    elif sys.argv[iarg].upper() == '-MINEL':  # if min elevation 
        iarg = iarg + 1
        lowel = float( sys.argv[iarg])
        print( "Using elevations > %7.2f (d) for Cold load calculation" % (lowel))
    elif sys.argv[iarg].upper() == '-D':
        print( 'Adding Debug Printing')
        doDebug = True
    else:
        break
    iarg = iarg + 1
    timearg = iarg
    namearg = iarg+1
# end of while not reading file names

#first argument is the averaging time in seconds
avetimesec = float(sys.argv[timearg])
print( "Average time: ", avetimesec, " (seconds)")
newObs = False
allFiles = False

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
maxvel = 180.
minvel = -550.
# currently used
maxvel = 200.
minvel = -180.

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

# now run through and find hot and cold loads obs
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
        print( 'File is not an astronomy file: ',filename)
        continue
    else:
        extension = parts[nparts-1]
    extension = extension.upper()
    if (extension != 'HOT') and (extension != 'AST') and (extension != 'CLD'):
        print( 'Extension not recognized : ', parts[nparts-1])
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
    print( "Found %3d Hot load observations" % nhot)
else:
    print( "No Hot load data, can not calibrate")
    exit()

xv = hot.xdata * 1.E-6
yv = hot.ydataA
nData = len( xv)        # get data length and indicies for middle of spectra 
n6 = int(nData/6)
n26 = int(2*n6)
n46 = int(4*n6)
n56 = int(5*n6)

xa = 200
#xb = 1024-xa
xb = nData-200


#previously just copy
hv = copy.deepcopy(yv)
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
    print( 'Min Vel at channel: ',xa, minvel)
    print( 'Max Vel at channel: ',xb, maxvel)
                                   
# will need smoothed hot load values in remaining calc
#for iii in range(1,(nData-2)):
#    hv[iii] = (yv[iii-2]+yv[iii-1]+yv[iii]+yv[iii+1]+yv[iii+2])/5.
#hv = yv

if doDebug:
    print( 'Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % (minGlat, maxGlat))

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

    if rs.telel < lowel:  #if elevation too low for a cold load obs
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
    print( "No high galactic latitude data: can not calibrate")
    exit()
else:
    high.ydataA = scalefactor * high.ydataA/nhigh
    print( "Found %3d High Galactic Latitude spectra" % (nhigh))
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
# if this was a new obs; restart the sums

def average_spec( ave_spec, in_spec, nave, firstutc, lastutc):
    """
    Averages spectra and deal with issues of weighting functions by observing duration
    """

    # if restarting the sum
    if nave == 0:
        ave_spec = copy.deepcopy(in_spec)  # initial spectrum is one just read
        firstutc = in_spec.utc
        lastutc = in_spec.utc
        nave = 1
        # print( 'Xmin: ', min(ave_spec.xdata)/1e6, 'Xmax: ', max(ave_spec.xdata),' MHz')
        # sums are weighted by durations
        ave_spec.ydataA = in_spec.ydataA * in_spec.durationSec
        # keep track of observing time for weighted sum
        ave_spec.durationSec = in_spec.durationSec
    else: # else not enough time yet, average ave_spec data
        lastutc = in_spec.utc
        ave_spec.count = ave_spec.count + in_spec.count
        nave = nave + 1
        ave_spec.ydataA = ave_spec.ydataA + (in_spec.ydataA * in_spec.durationSec)
        # keep track of observing time for weighted sum
        ave_spec.durationSec = ave_spec.durationSec + in_spec.durationSec
    return ave_spec, in_spec, nave, firstutc, lastutc

def compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold):
    nData = len(xv)
    vel = np.zeros(nData)
    tsys = np.zeros(nData)    # initialize arrays
    tsky  = np.zeros(nData)    # initialize arrays

    for jjj in range (0, nData):
        vel[jjj] = c * (nuh1 - xv[jjj])/nuh1
        tsys[jjj] = yv[jjj]/gain[jjj]
        tsky[jjj] = tsys[jjj]
#        tsky[jjj] = tsys[jjj] - Tsys

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
        print( 'Imin Error computing baseline: ', imin)
        imin = 0
    if imin >= nData:
        print( 'Imin Error computing baseline: ', imin)
        imin = nData-1

    if imax < 0:
        print( 'Imax Error computing baseline: ', imax)
        imax = 0
    if imax >= nData:
        print( 'Imax Error computing baseline: ', imax)
        imax = nData-1

    ya = np.median(tsky[(xa-10):(xa+10)])
    yb = np.median(tsky[(xb-10):(xb+10)])
    slope = (yb-ya)/(xb-xa)
# This is the only difference between M and T processing
#    for iii in range( nData):
#        tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))

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

# compute median raw values
hotmedian = np.median(hv[n6:n56])
coldmedian = np.median(cv[n6:n56])
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

Tsys = np.median(trx[n6:n56])
tStdA = np.std(trx[n6:n26])
tStdB = np.std(trx[n46:n56])
print( "Median Receiver Temp: %7.2f +/- %5.2f (%5.2f %5.2f) (K)" % ( Tsys, (tStdA+tStdB)/2., tStdA, tStdB))

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
        print( "%3d: %s " % (nRead, filename))

    if nRead == nFiles:   # if this is the last file, must force output
        allFiles = True

    rs = radioastronomy.Spectrum()
#  print( filename)
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
# if not a sky observation
    if rs.telel < 0. and (not allFiles):
        continue

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
    lastdate = date

# if first time reading data, set obs parameters
    if lastfreq == 0.:
        lastfreq = rs.centerFreqHz 
        lastbw = rs.bandwidthHz
        lastgain = rs.gains[0]
        lastaz = rs.telaz
        lastel = rs.telel
        cold = copy.deepcopy( rs)
        firstutc = rs.utc
        lastutc = rs.utc
        ncold = 0

    if ncold > 0:
        # time difference is between mid-points of integrations
        dt = rs.utc - cold.utc 
        # add the time since midpoint of latests
        dt = dt + datetime.timedelta(seconds=rs.durationSec)
        lastdate = date

        newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
        newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl
        if newObs:
            if lastfreq != rs.centerFreqHz:
                print( "Change: LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz")
                lastfreq = rs.centerFreqHz
            if lastbw != rs.bandwidthHz:
                print( "Change: LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz")
                lastbw = rs.bandwidthHz
            if lastgain != rs.gains[0]:
                print( "Change: LastGain: ", lastgain, "New: ", rs.gains[0], " dB")
                lastgain = rs.gains[0]
            if newAzEl:
                print( "Change: LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees")
                lastaz = rs.telaz
                lastel = rs.telel

        # if this is the last file and there was not a change in observing setup 
        if allFiles and (not newObs):
            cold, rs, ncold, firstutc, lastutc = average_spec( cold, rs, ncold, firstutc, lastutc)

        # if time to average (or end of all files)
        if (dt > avetime) or newObs or allFiles:
            cold.ydataA = cold.ydataA/float(cold.durationSec)
            if doDebug:
                print( "Average duration: %7.1f " % (cold.durationSec))

            # not calibrating hot load observations.
            if cold.telel < 0.:
                # Reset the ncold count and restart sum
                ncold = 0
                continue

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

            tsky, vel = compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold)
            tSys = np.median(tsky[n6:n56])
            tStdA = np.std(tsky[n6:n26])
            tStdB = np.std(tsky[n46:n56])
            tStd = (tStdA+tStdB)/.2
            
            ymin = min(tsky[n6:n56])
            ymax = max(tsky[n6:n56])
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            # compute average time from first and last utcs
            aveutc,duration = radioastronomy.aveutcs( firstutc, lastutc)
            if doDebug and (efirstutc == lastutc):
                print( "first and last utcs are the same: ", firstutc)
            # keep the average time, but use the duration from the integration times.
            cold.utc = aveutc 
            # pull out coordinates for labeling
            az = cold.telaz
            el = cold.telel
            cold.azel2radec()    # compute ra,dec from az,el and average utc
            gallon = cold.gallon
            gallat = cold.gallat
            label = 'L,L=%5.1f,%5.1f' % (gallon, gallat)
            
            if doDebug:
                print( "First utc: ", firstutc)
                print( "Last  utc: ", lastutc)
                print( "Ave   utc: ", aveutc)
            avedatetime = cold.utc.isoformat()
            datestr = avedatetime.split('T')
            atime = datestr[1]
            timeparts = atime.split('.')
            labeltime = timeparts[0]
            label = '%s, A,E: %5s,%5s, L,L: %5.1f,%5.1f' % (labeltime, az, el, gallon, gallat)
            if minel == maxel:
                label = '%s L,L=%5.1f,%5.1f (%d)' % (labeltime, gallon, gallat, ncold)
            else:
                label = '%s L,L=%5.1f,%4.1f A,E=%4.0f,%4.0f' % (labeltime, gallon, gallat, az, el)
            print( ' Max: %9.1f  Median: %8.2f +/- %5.2f %3d %s' % (ymax, tSys, tStd, ncold, label))
            if gallat < 7.5 and gallat > -7.5:
                plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
            elif gallat < 15. and gallat > -15.:
                plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
            else:
                plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            nplot = nplot + 1
            ncold = 0

    # end if a cold file
    if allFiles:
        break

    if rs.telel > 0:
    # Average in most recently read spectrum
        cold, rs, ncold, firstutc, lastutc = average_spec( cold, rs, ncold, firstutc, lastutc)

    #end for all files to sum
# end of all files to read

if doDebug:
    print( 'Number of remaining observations not plotted: ', ncold)

# if observations cover multiple days
if firstdate != lastdate:
    date = firstdate + " - " + lastdate
else:
    date = firstdate

if myTitle == "":
    myTitle = "%s    " % (date)
else:
    myTitle = myTitle + "  "

if (firstaz == otheraz):
    myTitle = myTitle + "Az = %6.1f, " % (firstaz)
else:
    myTitle = myTitle + "Az = %6.1f to %6.1f, " % (firstaz, otheraz)

if minel == maxel:
    myTitle = myTitle + " El=%6.1f" % (minel)
else:
    myTitle = myTitle + " El=%6.1f to %6.1f" % (minel, maxel)

fig.canvas.set_window_title(myTitle)
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
#plt.xlim(-600., 300.)
#plt.xlim(-350., 350.)
plt.xlim(-250., 250.)
# keep the plot from becoming too narrow
dy = yallmax - yallmin
if dy < 8:
    dy = 8
# set the y scale to go above and below all data
plt.ylim((yallmin-(dy/8.)), (yallmax+(dy/4.)))
plt.title(myTitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.show()
