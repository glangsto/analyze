1#Python Script to plot calibrated/baseline-fit  NSF record data.
#plot the raw data from the observation
#HISTORY
#20Apr16 GIL use external baseline fitting function
#20Apr09 GIL updated polynomial fit inputs, was computing 3rd order, not 2nd order
#20Apr06 GIL incorperate gain factor measurements as test before mapping
#20Apr04 GIL allow lower Galactic Latitude
#20Feb15 GIL normalize for number of spectra averaged
#20JAN03 GIL compute 20 RMSs across the spectrum and take the median of these as representative
#20JAN02 GIL compute RMS in signal free region; test for gain correction
#19DEC30 GIL add title option
#19OCT08 GIL update constant values and comments
#19OCT03 GIL add options to fit a polynomial baseline and setting velocity ranges
#19OCT03 GIL add polynomical baseline fit option
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
import datetime as dt
import radioastronomy
import copy
import interpolate
import gainfactor as gf

try:
    from PyAstronomy import pyasl
    baryCenterAvailable = True
except:
    print("!!!! Unable to import PyAstronomy !!!!")
    print("Can not compute Bary Center velocity offset")
    print("In Linux, try: ")
    print("sudo pip install PyAstronomy")
    baryCenterAvailable = False

try:
    import jdutil
except:
    print("!!!! Unable to import jdutil !!!!")
    print("This python file (jdutil.py) must be in your python path")
    print("Can not compute Julian Date from datetime")
    baryCenterAvailable = False

# default values
avetimesec = 3600.
# put your list of known RFI features here.  Must have at least two.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [7, 7]

nargs = len(sys.argv)
#first argument is the averaging time in seconds
timearg = 1
namearg = 2

if nargs < 3:
    print("M: Median baseline calibrated horn observations")
    print("Usage: M [-F] [-L <velocity>] [-H <velocity>] <average_seconds> <files>")
    print("Where <average_seconds>: Number of seconds of observations to average.")
    print("-F optionally do a polynomial baseline fit")
    print("-L optionally set the low velocity region for baseline fit")
    print("-H optionall set the high velocity region for baseline fit")
    print("Observation file list must include at least one hot load file")
    print("")
    print("Glen Langston - NSF   November 22, 2019")
    exit()

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
# default baseline fitting is 2nd order
fitOrder = 2
# optionally turn on debug plotting
doPoly = False
doRa = True
doDebug = False
nGain = 20
# processor index for gain update; zero means no gain update
pIndex = 0
firstRun = True
myTitle = ""      # Default no input title
maxvel = 220.
minvel = -maxvel

iarg = 1
# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-S':
        print('Folding specectra')
        doFold = True
    elif sys.argv[iarg].upper() == '-F':
        iarg = iarg+1
        fitOrder = int( sys.argv[iarg])
        if fitOrder < 0:
            fitOrder = 0
        elif fitOrder > 10:
            fitOrder = 10
        doPoly = True
        if fitOrder == 0:
            print("Fitting a constant baseline")
        elif fitOrder == 1:
            print("Fitting a linear baseline")
        else:
            print("Fitting a %d-nd order polynomical baseline" % (fitOrder))
    elif sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-G':
        iarg = iarg+1
        nGain = int( sys.argv[iarg])
        print('Dividing the Spectrum into %3d chunks and computing RMS for each' % (nGain))
    elif sys.argv[iarg].upper() == '-I':
        iarg = iarg+1
        pIndex = int( sys.argv[iarg])
        print('Updating Gain using processor Index %d' % (pIndex))
    elif sys.argv[iarg].upper() == '-L':
        iarg = iarg+1
        minvel = np.float( sys.argv[iarg])
        print('Minium (low)  velocity for sum: %7.2f km/sec' % (minvel))
    elif sys.argv[iarg].upper() == '-H':
        iarg = iarg+1
        maxvel = np.float( sys.argv[iarg])
        print('Maximum (high) velocity for sum: %7.2f km/sec' % (maxvel))
    elif sys.argv[iarg].upper() == '-MINEL':  # if min elevation 
        iarg = iarg + 1
        lowel = float( sys.argv[iarg])
        print("Using elevations > %7.2f (d) for Cold load calculation" % (lowel))
    elif sys.argv[iarg].upper() == '-D':
        print('Adding Debug Printing')
        doDebug = True
    elif sys.argv[iarg].upper() == '-T':   # if plot title provided
        iarg = iarg+1
        myTitle = sys.argv[iarg]
        print('Plot Title : ', myTitle)
    else:
        break
    iarg = iarg + 1
    timearg = iarg
    namearg = iarg+1
# end of while not reading file names

#first argument is the averaging time in seconds
try:
    avetimesec = float(sys.argv[timearg])
except:
    print("Error parsing time: %s" % (sys.argv[timearg]))
    avetimesec = 3600.
    # assume time is actually a file name
    namearg = timearg

print("Average time: ", avetimesec, " (seconds)")
newObs = False
allFiles = False

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-c','-m','-y','-c','-m','-y','-c','-m','-y','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
xallmax = -9.e9
xallmin = 9.e9
yallmax = -9.e9
yallmin = 9.e9

c = 299792.458  # (v km/sec)
nuh1 = 1420.4056E6 # neutral hydrogen frequency (Hz)
thot = 285.0  # define hot and cold
tcold = 10.0
tmin = 20.0 
tmax = 999.0 # define reasoanable value limits

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
        # print 'Xmin: ', min(ave_spec.xdata)/1e6, 'Xmax: ', max(ave_spec.xdata),' MHz'
        # sums are weighted by durations
        ave_spec.ydataA = (in_spec.ydataA/in_spec.count) * in_spec.durationSec
        # keep track of observing time for weighted sum
        ave_spec.durationSec = in_spec.durationSec
    else: # else not enough time yet, average ave_spec data
        lastutc = in_spec.utc
        ave_spec.count = ave_spec.count + in_spec.count
        nave = nave + 1
        ave_spec.ydataA = ave_spec.ydataA + (in_spec.ydataA * in_spec.durationSec/in_spec.count)
        # keep track of observing time for weighted sum
        ave_spec.durationSec = ave_spec.durationSec + in_spec.durationSec
    return ave_spec, in_spec, nave, firstutc, lastutc

def compute_vbarycenter( spectrum):
    """ 
    Compute the velocity correction to Barycentric for this date and direction
    """
    global firstRun 

    if baryCenterAvailable: 
        longitude = spectrum.tellon
        latitude = spectrum.tellat
        altitude = spectrum.telelev
        ra2000 = spectrum.ra
        dec2000 = spectrum.dec

        # need to convert date time to Julian date
        jd = jdutil.datetime_to_jd(spectrum.utc)

# Calculate barycentric correction (debug=True show
# various intermediate results)
        corr, hjd = pyasl.helcorr(longitude, latitude, altitude, \
                                      ra2000, dec2000, jd, debug=doDebug)
        if doDebug or firstRun:
            print("Barycentric correction [km/s]: %8.3f" % (corr))
            firstRun = False
    else:
        corr = 0.
    return corr

def velocity_to_indicies( vel, minvel, maxvel):
    """
    Function to compute indices from velocity array and target velocities
    """

    nData = len(vel)
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
        print('Imin Error computing baseline: ', imin)
        imin = 0
    if imin >= nData:
        print('Imin Error computing baseline: ', imin)
        imin = nData-1

    if imax < 0:
        print('Imax Error computing baseline: ', imax)
        imax = 0
    if imax >= nData:
        print('Imax Error computing baseline: ', imax)
        imax = nData-1

# some channels will be selected using x[imin:imax] so,
# now make sure range increases
    if imin > imax:
        temp = imin
        imin = imax
        imax = temp
    return imin, imax
    
def compute_tsky_hotcold( yv, gain, vel, minvel, maxvel):
    """
    Compute the System Temperature (Tsky) from the hot and cold load observations
    where:
    spectrum     Spectrum of interest in Counts
    gain         Precalibrated gain array for scaling the Spectral counts
    minvel       Minimum Velocity to sum
    """
    nData = len(vel)
    tsys = np.zeros(nData)    # initialize arrays

    if nData < 16:
        print("Error, velocity array too small: ",nData)

    for jjj in range(nData):
        tsys[jjj] = yv[jjj]/gain[jjj]

    if flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
        icenter = int(nData/2)
        tsys[icenter] = (tsys[icenter-3] + tsys[icenter+3])*.5
        tsys[icenter-1] = (3.*tsys[icenter-3] + tsys[icenter+3])*.25
        tsys[icenter+1] = (tsys[icenter-3] + 3.*tsys[icenter+3])*.25
        tsys[icenter-2] = (2.*tsys[icenter-3] + tsys[icenter+3])*.3333
        tsys[icenter+2] = (tsys[icenter-3] + 2.*tsys[icenter+3])*.3333

    # get indicies for given max and min velocities
    imin, imax = velocity_to_indicies( vel, minvel, maxvel)

    nfit = 10
    # if fitting a polynomial
    if doPoly:
        baseline = gf.fit_baseline( vel, tsys, imin, imax, 2*nfit, fitOrder)
        tsys = tsys - baseline
    #        tsys = baseline  ; as a test, plotted only the baseline fit
    else: # else a linear baseline
        ya = np.median(tsys[(imin-nfit):(imin+nfit)])
        yb = np.median(tsys[(imax-nfit):(imax+nfit)])
        slope = (yb-ya)/(imax-imin)
#                baseline = tsys
#                print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
# This is the only difference between M and T processing
        if doDebug:
            print("Slope %7.3f K/channel" % (slope))
        for iii in range( nData):
            tsys[iii] = tsys[iii] - (ya + (slope*(iii-imin)))
        
    return tsys

# initialize counters for processing many spectra
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
DT = dt.timedelta(seconds=0.)

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
        print('File is not an astronomy file: ',filename)
        continue
    else:
        extension = parts[nparts-1]
    extension = extension.upper()
    if (extension != 'HOT') and (extension != 'AST') and (extension != 'CLD'):
        print('Extension not recognized : ', parts[nparts-1])
        continue

    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

    if rs.telel < 0:
        if nhot == 0:
            hot = copy.deepcopy( rs)
            hot.ydataA = rs.ydataA/rs.count
            nhot = 1
        else:
            hot.ydataA = hot.ydataA + (rs.ydataA/rs.count)
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
    hot.ydataA = hot.ydataA / float(nhot)
    print("Found %3d Hot load observations" % nhot)
else:
    print("No Hot load data, can not calibrate")
    exit()

#xv = hot.xdata * 1.E-6
xv = hot.xdata
yv = hot.ydataA
nData = len( xv)        # get data length and indicies for middle of spectra 
n6 = int(nData/6)
n56 = 5*n6
rmss = np.zeros( nGain)
dx = int( nData/(nGain+4))
print('Computing %d RMS in %d channels each' % (nGain, dx))

#previously just copy
hv = copy.deepcopy(yv)
if flagRfi:
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

vel = np.zeros(nData)
# create index array
for jjj in range (nData):
    vel[jjj] = c * (nuh1 - xv[jjj])/nuh1

xa, xb = velocity_to_indicies( vel, minvel, maxvel)
print('Min Vel, channel: ',xa, minvel)
print('Max Vel, channel: ',xb, maxvel)

#if doDebug:
if True:
    print('Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % (minGlat, maxGlat))

# assume only a limited range of galactic latitudes are available
# not range about +/-60.
use60Range = False

# all galactic latitudes above minimum can be used for calibration
if (minGlat < -50.) or (maxGlat > 50.):
    minGlat = -50.
    maxGlat = 50.
else: # else no high galactic latitude data
    # use highest galactic latitudes - +/-5.degrees
    if -minGlat > maxGlat:  # if negative latitudes higher
        minGlat = minGlat + 5.
        maxGlat = 90.
    else: # else positive latitudes higher
        maxGlat = maxGlat - 5.
        minGlat = -90.

if True:
    print('Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % (minGlat, maxGlat))

# but if no low latitude data are available use minimum
# now average coldest data for calibration
rs = radioastronomy.Spectrum()
for filename in names:

    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

    if rs.telel < lowel:  #if elevation too low for a cold load obs
        continue

    if (rs.gallat > maxGlat) or (rs.gallat < minGlat):
        if doDebug:
            print ("%s: %7.1f %7.1f" % (filename, rs.gallon, rs.gallat))
        if nhigh == 0:
            high = copy.deepcopy( rs)
            high.ydataA = rs.ydataA/rs.count
            nhigh = 1
        else:
            high.ydataA = high.ydataA + (rs.ydataA/rs.count)
            high.count = high.count + rs.count
            nhigh = nhigh + 1

if nhigh < 1.:
    print("No high galactic latitude data: can not calibrate")
    exit()
else:
    high.ydataA = high.ydataA/nhigh
    print("Found %3d High Galactic Latitude spectra" % (nhigh))
    yv = high.ydataA

    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

# finally compute gain on a channel by channel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)
#    gain[iii] = hv[iii]/thot
    vel[iii] = c * (nuh1 - xv[iii])/nuh1

#now want gain using only hot counts, but must add in Tsys
tSys = cv/gain
tSysMiddle = np.median(tSys[n6:n56])
CountsPerKelvin = np.median(gain[n6:n56])

# for remainder of calculations only use hot counts for calibration
for iii in range(nData):
    gain[iii] = hv[iii]/(thot + tSysMiddle - tcold)

fig, ax1 = plt.subplots(figsize=(10, 6))

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

# compute the reciever temperature
trx = np.zeros(nData)
zeros = np.zeros(nData)
for iii in range(nData):
    trx[iii] = (cv[iii]/gain[iii]) - tcold

Tsys = np.median(trx[xa:xb])
print("Median Receiver Temp: %7.2f (K)" % ( Tsys))

avetime = dt.timedelta(seconds=avetimesec)
nRead = 0          # so far no file names read

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
        print("%3d: %s " % (nRead, filename))

    if nRead == nFiles:   # if this is the last file, must force output
        allFiles = True

    rs = radioastronomy.Spectrum()
#  print filename
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
        veryfirstutc = rs.utc
        verylastutc = rs.utc
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

    DT = dt.timedelta(seconds=rs.durationSec)

    newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
    newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl

    # time difference is between mid-points of integrations
    DT = DT + (rs.utc - cold.utc)
    # add the time since midpoint of latests
    lastdate = date
    if rs.utc > verylastutc:
        verylastutc = rs.utc
    if rs.utc < veryfirstutc:
        veryfirstutc = rs.utc

    # if this is the last file and there was not a change in observing setup 
    if allFiles and (not newObs):
        cold, rs, ncold, firstutc, lastutc = average_spec( cold, rs, ncold, firstutc, lastutc)

    # if time to average (or end of all files)
    if (DT > avetime) or newObs or allFiles:
        cold.ydataA = cold.ydataA/float(cold.durationSec)
        if doDebug:
            print("Average duration: %7.1f " % (cold.durationSec))

        if cold.telel < 0.:
            ncold = 0
            continue

        yv = cold.ydataA
        if flagRfi:
            yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
        xmin = min(xv[xa:xb])
        xmax = max(xv[xa:xb])
        xallmin = min(xmin, xallmin)
        xallmax = max(xmax, xallmax)
        count = cold.count
        note = cold.noteA

        ncolor = min(nmax-1, nplot)  # select color for this spectrum

        # compute velocity correction for this direction and date
        corr = gf.compute_vbarycenter( cold)
        velcorr = vel + corr

        tsky = compute_tsky_hotcold( yv, gain, velcorr, minvel, maxvel)
        ymed = np.median(tsky[xa:xb])
        ystd = np.std(tsky[xa:xb])

        ia = 2*dx
        ib = 3*dx
        # divide spectrum into bins and compute RMS for each
        for lll in range(nGain):
            rmss[lll] = np.std(tsky[ia:ib])
            ia=ia+dx
            ib=ib+dx

        rmsmin = min( rmss)
        rmsmax = max( rmss)
        rmsmed = np.median( rmss)

        # compute average time from first and last utcs
        aveutc,duration = radioastronomy.aveutcs( firstutc, lastutc)
        if doDebug and (efirstutc == lastutc):
            print("first and last utcs are the same: ", firstutc)
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
            print("First utc: ", firstutc)
            print("Last  utc: ", lastutc)
            print("Ave   utc: ", aveutc)
        # create label string based on the average time
        avedatetime = cold.utc.isoformat()
        datestr = avedatetime.split('T')
        atime = datestr[1]
        timeparts = atime.split('.')
        labeltime = timeparts[0]
        label = '%s, A,E: %5.1f,%5.1f, L,L: %4.0f,%4.0f' % (labeltime, az, el, gallon, gallat)
        if minel == maxel:
            label = '%s L,L=%5.1f,%5.1f (%d)' % (labeltime, gallon, gallat, ncold)
        else:
            label = '%s L,L=%5.1f,%5.1f A,E=%3.0f,%3.0f' % (labeltime, gallon, gallat, az, el)
        if doRa:
            label = label + ' R,D=%5.1f,%5.1f' % (cold.ra, cold.dec)

        if pIndex > 0:
            gainfactor = gf.compute_gain_factor( pIndex, aveutc, el)
            print("Gain Factor: %7.2f; for Index %d, %s, %7.1f" % (gainfactor, pIndex, aveutc, el))
            tsky = gainfactor * tsky

        ymin = min(tsky[xa:xb])
        ymax = max(tsky[xa:xb])
        yallmin = min(ymin,yallmin)
        yallmax = max(ymax,yallmax)

        if nplot < 1:
            print('    Max     Median              RMSs           N      Time   Galactic ')
            print('    (K)      (K)      Median    Bottom    Top  Ave            Lon, Lat ')
        print('%9.1f %8.2f %8.2f %8.2f %8.2f %3d %s' % (ymax, ymed, rmsmed, rmsmin, rmsmax, ncold, label))
        if gallat < 7.5 and gallat > -7.5:
            lw = 4
        elif gallat < 15. and gallat > -15.:
            lw = 2
        else:
            lw = 1
        plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=lw)
        nplot = nplot + 1
        ncold = 0

    if newObs:
        if lastfreq != rs.centerFreqHz:
            print("Change: LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz")
            lastfreq = rs.centerFreqHz
        if lastbw != rs.bandwidthHz:
            print("Change: LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz")
            lastbw = rs.bandwidthHz
        if lastgain != rs.gains[0]:
            print("Change: LastGain: ", lastgain, "New: ", rs.gains[0], " dB")
            lastgain = rs.gains[0]
        if newAzEl:
            print("Change: LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees")
            lastaz = rs.telaz
            lastel = rs.telel
    # end if a cold file
    if allFiles:
        break

    if rs.telel > 0:
    # Average in most recently read spectrum
        cold, rs, ncold, firstutc, lastutc = average_spec( cold, rs, ncold, firstutc, lastutc)

    if ncold == 1:
        # At each sum restart, must recompute velocity axis
        xv = cold.xdata
        vel = copy.deepcopy( xv)
        vel = nuh1 - vel         # converting to frequency offset
        vel = (c/nuh1) * vel     # now conver to velocity
        lastutc = cold.utc

    #end for all files to sum
# end of all files to read

if doDebug:
    print('Number of remaining observations not plotted: ', ncold)

if veryfirstutc > verylastutc:
    temp = veryfirstutc
    veryfirstutc = verylastutc
    temp = verylastutc

utcstr = str(veryfirstutc)
parts = utcstr.split(' ')
date = parts[0]
nd = len(date)
firstdate = date[2:nd]

utcstr = str(verylastutc)
parts = utcstr.split(' ')
date = parts[0]
nd = len(date)
lastdate = date[2:nd]

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
plt.xlim(minvel, maxvel)
# keep the plot from becoming too narrow
if yallmax < 8:
    yallmax = 8
# set the y scale to go above and below all data
plt.ylim((yallmin*.97)-1., 1.15*yallmax)
plt.title(myTitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.show()
