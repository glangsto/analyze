#Python Script to plot calibrated/baseline-fit  NSF record data.
#plot the raw data from the observation
#HISTORY
#19OCT03 GIL add options to fit a polynomial baseline and setting velocity ranges
#19SEP26 GIL merge m.py and sumhotcold functions
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime as dt
import radioastronomy
import copy
import interpolate
import Moment
import os
try:
    #            from __future__ import print_function, division
    from PyAstronomy import pyasl
    baryCenterAvailable = True
except:
    print "!!!! Unable to import PyAstronomy !!!!"
    print "Can not compute Bary Center velocity offset"
    print "In Linux, try: "
    print "sudo pip install PyAstronomy"
    baryCenterAvailable = False

try:
    import jdutil
except:
    print "!!!! Unable to import jdutil !!!!"
    print "Can not compute Julian Date from datetime"
    baryCenterAvailable = False

# default values
avetimesec = 3600.
# put your list of known RFI features here.  Must have at least two.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]

nargs = len(sys.argv)
#first argument is the averaging time in seconds
timearg = 1
namearg = 2

# special printout on first execution
firstRun = True

# some special printout on first run
firstrun = True
# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False
# put list of RFI features here, for interpolation later
flagRfi = False
flagRfi = True
# to address an old problem, optionally allow folding spectra
doFold = False
# specify lowest  elevation for cold load averge
lowel = 60.
# optionally turn on debug plotting
doDebug = False
doPlot = False
doPoly = False
iarg = 1
maxvel = 200.
minvel = -maxvel

if nargs < 2:
    print "SUMHOTCOLD: Summarize a set of observations, computing the intensity integrals"
    print "usage: SUMHOTCOLD [-D] <aveTimeSeconds> <hot file> <cold file> <file names>"
    print "where"
    print " -D                  optionally print debug info"
    print " -P                  optionally plot spectra"
    print " -L <velocity>       optionally set the low velocity range for integration"
    print " -H <velocity>       optionally set the high velocity range for integration"
    print "where hot  file      name of the calibration hot  file for this observation"
    print "where cold file      name of the calibratino cold file for this observation"
    print "where aveTimeSeconds time range of data to be averaged before output, in seconds"
    print "where file names     list of files to examine. Should not include hot files"
    print ""
    exit()

# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-D':
        print 'Adding Debug Printing'
        doDebug = True
    elif sys.argv[iarg].upper() == '-F':
        print 'Fitting Parabolic Baseline'
        doPoly = True
    elif sys.argv[iarg].upper() == '-P':
        doPlot = True
    elif sys.argv[iarg].upper() == '-L':
        iarg = iarg+1
        minvel = np.float( sys.argv[iarg])
        print 'Minium (low)  velocity for sum: %7.2f km/sec' % (minvel)
    elif sys.argv[iarg].upper() == '-H':
        iarg = iarg+1
        maxvel = np.float( sys.argv[iarg])
        print 'Maximum (high) velocity for sum: %7.2f km/sec' % (maxvel)
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
    print "Error parsing time: %s" % (sys.argv[timearg])
    avetimesec = 3600.
    # assume time is actually a file name
    namearg = timearg

# hot load file and cold load files are first two file names
hotfilename = sys.argv[namearg]
coldfilename = sys.argv[namearg+1]
# now skip first two files
firstfilearg = namearg+2
names = sys.argv[firstfilearg:]
names = sorted(names)
nFiles = len(names)

avetime = dt.timedelta(seconds=avetimesec)
print "Average time: ", avetimesec, " (seconds)"

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
eloffset = -2.75
azoffset = -34.25
# latest offset after averaging previous best with current best
eloffset = 0.0
azoffset = 0.0

c = 299792.  # (v km/sec)
nuh1 = 1420.4056E6 # neutral hydrogen frequency (Hz)
thot = 285.0  # define hot and cold
tcold = 10.0
tmin = 40.0 
Tmax = 999.0 # define reasoanable value limits

# prepare to remove a linear baseline
nplot = 0
nsum = 0
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
DT = dt.timedelta(seconds=0.)

# Read the hot and cold files
hot = radioastronomy.Spectrum()
hot.read_spec_ast(hotfilename)
xv = hot.xdata
yv = copy.deepcopy( hot.ydataA)

nData = len( xv)        # get data length and indicies for middle of spectra 
if nData < 1:
    print "Error reading hot file: %s" %s (hotfilename)
    print "No data!"
    exit()

n6 = int(nData/6)
n56 = int(5*n6)

# default indexs for summing and interplating
xa = int(n6)
xb = int(n56)

#previously just copy
hv = copy.deepcopy(yv)
if flagRfi:
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

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

# now make sure range increases
    if imin > imax:
        temp = imin
        imin = imax
        imax = temp
    return imin, imax

vel = np.zeros(nData)
# need to assume all spectra have the same x axis
for jjj in range (nData):
    vel[jjj] = c * (nuh1 - xv[jjj])/nuh1

xa, xb = velocity_to_indicies( vel, minvel, maxvel)

hotMedian = np.median(hv[xa:xb])

if doDebug:
    print 'Min Vel at channel: ',xa, minvel
    print 'Max Vel at channel: ',xb, maxvel
                                   
cold = radioastronomy.Spectrum()
cold.read_spec_ast(coldfilename)
cv = copy.deepcopy(cold.ydataA)
nCold = len(cold.ydataA)

if nCold < 1:
    print "Error reading cold file: %s" % (coldfilename)
    print "No data!"
    exit()

if flagRfi:
    cv = interpolate.lines( linelist, linewidth, xv, cv) # interpolate rfi
coldMedian = np.median( cv[xa:xb])

# Compute gain on a channel by channel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)
    if gain[iii] <= 0.:
        print "Gain(%d): %7.2f %7.2f %7.2f %7.1f %7.1f" % (iii, gain[iii], hv[iii], cv[iii], thot, tcold)

#now want gain using only hot counts, but must add in Tsys
tSys = cv/gain
tSysMiddle = np.median(tSys[xa:xb])
CountsPerKelvin = np.median(gain[xa:xb])

# for remainder of calculations only use hot counts for calibration
for iii in range(nData):
    gain[iii] = hv[iii]/(thot + tSysMiddle - tcold)

def fit_baseline( xs, ys, imin, imax):

    xfit = np.concatenate( (xs[imin-20:imin+20],xs[imax-20:imax+20]))
    yfit = np.concatenate( (ys[imin-20:imin+20],ys[imax-20:imax+20]))
# calculate polynomial
    z = np.polyfit(xfit, yfit, 3)
    f = np.poly1d(z)

# calculate new x's and y's
    yout = f(xs)

    return yout

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
                                      ra2000, dec2000, jd, debug=False)
        if doDebug or firstRun:
            print "Barycentric correction [km/s]: %8.3f" % (corr)
            print "Heliocentric Julian day      : %9.2f" % (hjd)
            print "Lon, Lat                     : %8.3f,%8.3f" % (longitude, latitude)
            print "Ra, Dec (J2000)              : %8.3f,%8.3f" % (ra2000, dec2000)
            firstRun = False
    else:
        corr = 0.
    return corr

def write_sum( outname, xa, xb, ra, dec, gallon, gallat, time, vel, tsky, el, ratio, n):
    """ 
    write out the integral of an array over a relevant velocity range
    where
    outname   File to record sums
    xa        lower index of region to sum
    xb        upper index of region to sum
    ra,dec    Average Right Ascension, Declination coordinates of sum
    gallon,gallat    Average Galactic coordinates of sum
    time      String of average utc of sum
    vel       Array of velocities of sum
    tsky      Array of calibrated intensities in Kelvins
    el        Elevation of observation
    curMedian Median of current (raw) counts of this observation
    coldMedian  Median of the single cold load used to calibrate these observations
    n         Number of spectra averaged to produce this sum
    """
    global firstrun
    nvel = len(vel)
    nsky = len(tsky)
    if (nsky > nvel):
        print "Error writing sum: Velocites > Intensity array sizes: ", nvel, nsky

    if nvel < 5:
        print 'Not enough values in arrays to compute reliable sums'

    # compute reference velicity and increments
    nvel2 = int(nvel/2)
    dv = (vel[xa]-vel[xb])/np.float(xa-xb)

    if firstrun:
        currentdir = os.getcwd()
        #print 'Velocity resolution : %9.3f km/sec' % dv 
        print 'Current Directory: %s' % (currentdir) 

# get rms in region expected not to have signal
    vs = vel[xa:xa+10]
    ts = tsky[xa:xa+10]
    (ave,adev,rms1,var,skew,curt) = Moment.moment(vs, ts)    
    vs = vel[xb-10:xb]
    ts = tsky[xb-10:xb]
    (ave,adev,rms2,var,skew,curt) = Moment.moment(vs, ts)    
    rms = 0.5*(rms1+rms2)

# now get the sum over the range of interest
    vs = vel[xa:xb]
    ts = tsky[xa:xb]
    (ave,adev,sdev,var,skew,curt) = Moment.moment(vs, ts)    
    tmax = max(ts)
    tstd = sdev

# now compute the intensity weighted velocity sum.
    vts = vel[xa:xb]*tsky[xa:xb]

    # compute the weighted average velocity
    (vave,vadev,vsdev,vvar,vskew,vcurt) = Moment.moment(vs, vts)    
    # normalize for average temperature
    if ave == 0: 
        print 'Average temperature is zero, can not normalize average velocity'
    else:        
        vave = vave/ave
        vsdev = vave/ave

    if not os.path.isfile(outname):
        f = open( outname, 'w')
#        f.write( '#  LON     LAT      T_sum     T_std     V_ave     V_std    Time \n')
#       f.write( '#   RA     DEC     LON     LAT      T_sum     T_std     V_ave     V_std  N  Time \n')
        f.write( '#   RA     DEC     LON     LAT      T_sum +/- T_std   T_max  V_ave+/-V_std  El    N  Time \n')
#        f.write( '#   RA     DEC     LON     LAT      T_sum     T_max     T_std   V_ave  V_std  El    N  Time \n')
        # Now label the terminal outout
    else:
        f = open( outname, 'a+')

    if firstrun:
        print '   Galactic         Intensity         Intensity     Count    '
        print '   Lon, Lat       Sum      Rms       Peak   Rms     Ratio  Elev  Date'
        print '    (d, d)          (K km/sec)           (K)               (d)   '

    firstrun = False
# File = 2017-01-03.sum
# RA  DEC   T_sum  T_std   V_ave  V_std   Time  
# 54.56  -35.32    -5.587     3.039  -624.929   111.851 20:02:38
    tsum = 0
    # integer number of samples
    dx = xb-xa
    for iii in range(dx):
        tsum = tsum+ts[iii]
# take absolute value of velocity
    if dv < 0.:
        dv = -dv
    tsum = tsum * dv  # scale integral by velocity resolution
    sdev = adev * dv * np.sqrt(float(dx)) # scale error estimate by approriate number of channels
    if tsum < -3.*sdev:
        print 'ave: %7.2f dv=%7.2f  tsum=%7.2f ave*dv=%7.2f' % (ave, dv, tsum, ave*dv*float(dx))

    else:
        print '%6.1f,%5.1f %9.2f +/- %6.2f %6.1f +/- %5.1f %6.3f %4.0f %s' % \
            (gallon, gallat, tsum, sdev, tmax, tstd, ratio, el, time)
        vsdev = rms * dv * np.sqrt(float(dx)*abs(sdev/tsum)) # scale error estimate by approriate number of channels
        label = '%7.2f %7.2f %7.2f %7.2f %9.2f %9.2f %7.2f %7.2f %6.2f %4.1f %3d %s\n' % \
                (ra, dec, gallon, gallat, tsum, sdev, tmax, vave, vsdev, el, n, time)
        f.write(label)
    f.close()

def model_tsys_el( tSys90, el):
    """  
    Model the system temperature versus elevation for input elevation
    tSys90    System temperature at Zenith
    el        elevation of telescope degrees
    returns the estimated System temperature at this elevation
    """
    pi = 3.1415926
    toradians = pi/180.
    secz = 1.
    if el > 80.:
        secz = 1.
    elif el > 0.1:
        secz = 1./np.cos(toradians*el)

    tSky = tcold * secz
    if tSky > Tmax:
        print "Tsky: %9.1f el: %7.2f secz=%7.3f" % (tSky, el, secz)

    if (el < 45.) and (el > 0.):
        tGround = (thot/2.)*np.cos(el*toradians)**2.
    else:
        tGround = 0.
    tSysEl = np.float(tSys90 + tSky + tGround)
    ratio = np.float(tSysEl/tSys90)
    if doDebug:
        print "Tsys Estimate: el=%5.1f tSys90=%5.1f Tsky=%5.1f tGround=%5.1f ratio=%6.3f" % \
            (el, tSys90, tSky, tGround, ratio)
    return tSysEl
    
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
    return ave_spec, nave, firstutc, lastutc

def compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold):
    nData = len(xv)
    vel = np.zeros(nData)
    tsys = np.zeros(nData)    # initialize arrays
    tsky  = np.zeros(nData)    # initialize arrays

    for jjj in range (0, nData):
        vel[jjj] = c * (nuh1 - xv[jjj])/nuh1
        tsys[jjj] = yv[jjj]/gain[jjj]
        tsky[jjj] = tsys[jjj]

    if flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
        icenter = int(nData/2)
        tsky[icenter] = (tsky[icenter-2] + tsky[icenter+2])*.5
        tsky[icenter-1] = (3.*tsky[icenter-2] + tsky[icenter+2])*.25
        tsky[icenter+1] = (tsky[icenter-2] + 3.*tsky[icenter+2])*.25

    xa, xb = velocity_to_indicies( vel, minvel, maxvel)

    if doPoly:
        baseline = fit_baseline( xv, tsky, xa, xb)
        tsky = tsky - baseline
    else:
        ya = np.median(tsky[(xa-10):(xa+10)])
        yb = np.median(tsky[(xb-10):(xb+10)])
        slope = (yb-ya)/(xb-xa)
#                baseline = tsky
#                print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
# This is the only difference between M and T processing
        if doDebug:
            print "Slope %7.3f K/channel" % (slope)
        for iii in range( nData):
            tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))

    return tsky, vel

if doPlot:
    fig, ax1 = plt.subplots(figsize=(10, 6))
ymin = 10000000000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
nsum = 0
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
    lastdate = date

# if first time reading data, set obs parameters
    if lastfreq == 0.:
        lastfreq = rs.centerFreqHz 
        lastbw = rs.bandwidthHz
        lastgain = rs.gains[0]
        lastaz = rs.telaz
        lastel = rs.telel
        rsum = copy.deepcopy( rs)
        firstutc = rs.utc
        lastutc = rs.utc
        tSysEl = model_tsys_el( tSysMiddle, lastel)
        tSysFactor = tSysEl/tSysMiddle
        nsum = 0

    if nsum > 0:
        # time difference is between mid-points of integrations
        DT = rs.utc - rsum.utc 
        # add the time since midpoint of latests
        DT = DT + dt.timedelta(seconds=rs.durationSec)
        lastdate = date

        newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
        newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl
        if newObs:
            if lastfreq != rs.centerFreqHz:
                print "Change: LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz"
                lastfreq = rs.centerFreqHz
            if lastbw != rs.bandwidthHz:
                print "Change: LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz"
                lastbw = rs.bandwidthHz
            if lastgain != rs.gains[0]:
                print "Change: LastGain: ", lastgain, "New: ", rs.gains[0], " dB"
                lastgain = rs.gains[0]
            if newAzEl:
                print "Change: LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees"
                lastaz = rs.telaz
                lastel = rs.telel
            if minel > rs.telel:
                minel = rs.telel
            if maxel < rs.telel:
                maxel = rs.telel

        # if this is the last file and there was not a change in observing setup 
        if allFiles and (not newObs):
            rsum, nsum, firstutc, lastutc = average_spec( rsum, rs, nsum, firstutc, lastutc)

        # if time to average (or end of all files)
        if (DT > avetime) or newObs or allFiles:
            rsum.ydataA = rsum.ydataA/float(rsum.durationSec)
            if doDebug:
                print "Average duration: %7.1f " % (rsum.durationSec)

            if rsum.telel < 0.:
                nsum = 0
                continue

            xv = rsum.xdata
            yv = copy.deepcopy(rsum.ydataA)
            currentMedian = np.median( yv[xa:xb])
            ratio = (coldMedian/currentMedian) * tSysFactor
            if flagRfi:
                yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

            yv = yv*ratio # scale to constant system tmemperature
            xmin = min(xv)
            xmax = max(xv)
            xallmin = min(xmin, xallmin)
            xallmax = max(xmax, xallmax)
            count = rsum.count
            note = rsum.noteA
                    #print('%s' % note)
            ncolor = min(nmax-1, nplot) 


            tsky, vel = compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold)
            ymed = np.median(tsky[xa:xb])
            ystd = np.std(tsky[xa:xb])
            if ystd <= 0.:
                ystd = 0.001
            if ymed > Tmax:
                print "Tsys:%9.0f exceeds Tmax %7.2f" % (ymed, Tmax)
                nsum = 0
                continue
            if ystd > Tmax:
                print "Trms:%9.0f exceeds Tmax %7.2f" % (ystd, Tmax)
                nsum = 0
                continue

            ymin = min(tsky[xa:xb])
            ymax = max(tsky[xa:xb])
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            # compute average time from first and last utcs
            aveutc,duration = radioastronomy.aveutcs( firstutc, lastutc)
            if doDebug and (firstutc == lastutc):
                print "first and last utcs are the same: ", firstutc
            # keep the average time, but use the duration from the integration times.
            rsum.utc = aveutc 

            # compute the offset between earth's velocity relative to barycenter
            corr = compute_vbarycenter( rsum)
            velcorr = vel + corr

            # pull out coordinates for labeling
            az = rsum.telaz
            el = rsum.telel
            rsum.azel2radec()    # compute ra,dec from az,el and average utc
            ra = rsum.ra
            dec = rsum.dec
            gallon = rsum.gallon
            gallat = rsum.gallat
            label = 'L,L=%5.1f,%5.1f' % (gallon, gallat)
            
            if doDebug:
                print "First utc: ", firstutc
                print "Last  utc: ", lastutc
                print "Ave   utc: ", aveutc
            avedatetime = rsum.utc.isoformat()
            datestr = avedatetime.split('T')
            atime = datestr[1]
            timeparts = atime.split('.')
            labeltime = timeparts[0]
            label = '%s, A,E: %5s,%5s, L,L: %5.1f,%5.1f' % (labeltime, az, el, gallon, gallat)
            if minel == maxel:
                label = '%s L,L=%5.1f,%5.1f (%d)' % (labeltime, gallon, gallat, nsum)
            else:
                label = '%s L,L=%5.1f,%4.1f A,E=%4.0f,%4.0f' % (labeltime, gallon, gallat, az, el)
            if doDebug: 
                print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, (ymax-ymed)/ystd, nsum, label)

            hours = time[0:2]
            outname = date + '_' + hours + '.sum'
            # compute sums over a smaller range that the fit
            write_sum( outname, xa, xb, ra, dec, gallon, gallat, date + ' ' + time, \
                           velcorr, tsky, el, ratio, nsum)
            if (doPlot):
                lw = 1
                if gallat < 7.5 and gallat > -7.5:
                    lw = 4
                elif gallat < 15. and gallat > -15.:
                    lw = 2
                else:
                    lw = 1
                plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=lw)
                nplot = nplot + 1
            nsum = 0

        if nsum <= 0: # new observation, so recalculate the tSys Factor
            el = rs.telel
            tSysEl = model_tsys_el( tSysMiddle, el)
            tSysFactor = tSysEl/tSysMiddle
            nsum = 0

    # end if a cold file
    if allFiles:
        break

    if rs.telel > 0:
    # Average in most recently read spectrum
        rsum, nsum, firstutc, lastutc = average_spec( rsum, rs, nsum, firstutc, lastutc)

    if nsum == 1:
        # At each sum restart, must recompute velocity axis
        xv = rsum.xdata
        vel = copy.deepcopy( xv) # transfer frequency (Hz)
        vel = nuh1 - vel         # converting to frequency offset (Hz)
        vel = (c/nuh1) * vel     # now convert to velocity (km/sec)

    #end for all files to sum
# end of all files to read

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

if doPlot:
    fig.canvas.set_window_title(mytitle)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(14) 

    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14) 

if doPlot:
    plt.xlim(minvel, maxvel)
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
