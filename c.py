#Python Script to plot Calibrated NSF Horn Observations
#HISTORY
#19OCT08 GIL complete merger with M, T and SUMHOTCOLD
#19OCT03 GIL add options to fit a polynomial baseline and setting velocity ranges
#19OCT03 GIL add polynomical baseline fit option
#19SEP23 GIL use function for averaging 
#19SEP21 GIL fix exit on processing all files
#19SEP21 GIL fix exit on processing all files
#19SEP11 GIL do not write the Kelvins file until fixed (.kel)
#19JUN29 GIL debugginginfo added
#18DEC11 GIL add argument processing loop, saving products
#18DEC10 GIL initial version based on m.py
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime as dt
import radioastronomy
import hotcold
import copy
import interpolate
try:
    #            from __future__ import print_function, division
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

# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False # flag interpolate over spike in center of spectrum
doDebug = False   # flag printing debug info
#doDebug =True     # flag printing debug info
doSave = False    # flag saving intermediate files
flagRfi = True    # flag flagging RFI
doFold = False    # fold spectra to address an old issue; not normally used.
doSub = False     # define subtract baseline
doKelvins = False # save the kelvins calibration
doCalibrate = True # flagCalibrating observations
outDir = "./"     # define output directory for saving files
mytitle = ""      # define plot title
doHanning = False # hanning smooth hot load to reduce calibration noise

# put your list of known RFI features here.  Must have at least two, if flagRfi is true.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]   # integer number of channels to interpolate over
# currently used velocities for baseline fitting
maxvel = 300.
minvel = -300.
thot = 285.0  # define hot and cold in Kelvins
tcold = 10.0
nuh1 = 1420.40557177667 # neutral hydrogen frequency (MHz)
nureference = 1.E6*nuh1 # neutral hydrogen frequency (Hz)
# specify lowest  elevation for cold load averge
lowel = 60.
doPoly = False
doDebug = False
firstRun = True
maxvel = 200.
minvel = -maxvel
nFit = 10

nargs = len(sys.argv)
if nargs < 3:
    print("C: Calibrate Science Aficonado (NSF) horn observations")
    print("C: Optionally produce average hot and cold load spectra")
    print("Usage: C [options]  <average_seconds> <files>")
    print("")
    print("Where many parameters are optional:")
    print("-B Subtract a linear baseline fit to Spectra at Min and Max Velocities")
    print("   Min and max default velocities are: %7.1f, %7.1f km/sec" % (minvel, maxvel))
    print("-C Flag the Center channel, using interpolation.")
    print("   This removes a strong narrow feature created by many Software Defined Radios (SDRs)")
    print("-D Additional Debug printing.")
    print("-I Optionally Flag Known Radio Frequency Interference (RFI)")
    print("   Note you need to update the c.py program to add your list of Frequencies")
    print("   RFI frequencies are Location Dependent")
    print("-H <highvelocity>  High velocity for fit of baseline")
    print("-L <lowvelocity>   Low  velocity for fit of baseline")
    print("-HOT <hot load Temperature> Set the effective temperature of the hot load (Kelvins)")
    print("-N Not Calibrate.  This mode is used for tests of raw spectra")
    print("-O <output directory> Set the output directory for saved files")
    print("-R <Reference Frequency> Rest Frequency (Hz) used for Doppler calculations: %8.3f (MHz)" % (nureference*1.E-6))
    print("-S Save average spectra in files.  The Hot and Cold Load averages are saved, too.")
    print("   Average spectra have -ave added to their names")
    print("   Calibrated spectra have a .kel (for Kelvins) extension")
    print("-K Save the kelvins calibrated spectra as well")
    print("-X Hanning smooth the hot load observation to reduce calibration noise")
    print("-T <plot title String> Label for plot")
    print("-VA <low velocity> limit to exclude for baseline fitting")
    print("-VB <high velocity> limit to exclude for baseline fitting")
    print("Where:")
    print("   <average_seconds>: Number of seconds of observations to average.")
    print("   <average_seconds> is clock time, not observing time, so 3600. gives one plot for each hour")
    print("   <files> are Horn Observation files")
    print("   <files> must include both data pointed up (.ast) and down (.hot) observations")
    print("      All .hot files are assumed to have a system temperature of %7.1f K" % thot)
    print("")
    print(" -- Glen Langston (glangsto@nsf.gov), 2018 December 11")
    print("")
    exit()

#first argument is the averaging time in seconds
timearg = 1
namearg = 2
iarg = 1          # start searching for input flags
# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-F':
        doPoly = True
        print("Fitting a polynomial to the spectral baseline")
    elif sys.argv[iarg].upper() == '-I':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-MINEL':  # if min elevation 
        iarg = iarg + 1
        lowel = float( sys.argv[iarg])
        print("Using elevations > %7.2f (d) for Cold load calculation" % (lowel))
    elif sys.argv[iarg].upper() == '-N':   # if no calibration 
        print('Not calibrating Average Spectra in -ave files')
        doCalibrate = False
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-L':
        iarg = iarg+1
        minvel = np.float( sys.argv[iarg])
        print('Minium (low)  velocity for sum: %7.2f km/sec' % (minvel))
    elif sys.argv[iarg].upper() == '-H':
        iarg = iarg+1
        maxvel = np.float( sys.argv[iarg])
        print('Maximum (high) velocity for sum: %7.2f km/sec' % (maxvel))
    elif sys.argv[iarg].upper() == '-D':
        print('Adding Debug Printing')
        doDebug = True
    elif sys.argv[iarg].upper() == '-S':
        print('Saving Average Spectra in -ave files')
        doSave = True
    elif sys.argv[iarg].upper() == '-K':
        print('Saving Kelvins Calibration in .kel files')
        doSave = True
    elif sys.argv[iarg].upper() == '-B':
        print('Baseline subtraction')
        doSub = True
    elif sys.argv[iarg].upper() == '-X':
        print('Hanning Smooth the Hot load Observation')
        doHanning = True
    elif sys.argv[iarg].upper() == '-O':   # now look for flags with arguments
        iarg = iarg+1
        outDir = sys.argv[iarg]
        print('Output directory: ', outDir)
    elif sys.argv[iarg].upper() == '-HOT':   # now look for flags with arguments
        iarg = iarg+1
        thot = float(sys.argv[iarg])
        print('Hot Load Temperature: %7.2f K ', thot)
    elif sys.argv[iarg].upper() == '-R':   # reference frequency (Hz)
        iarg = iarg+1
        nureference = float(sys.argv[iarg])
        print('Reference Frequency: %8.3f (MHz): ', nureference * 1.E-6)
    elif sys.argv[iarg].upper() == '-T':   # now look for flags with arguments
        iarg = iarg+1
        mytitle = sys.argv[iarg]
        print('Plot title: ', mytitle)
    elif sys.argv[iarg].upper() == '-VA':   # now look for flags with arguments
        iarg = iarg+1
        minvel = float(sys.argv[iarg])
        print('Minimum velocity for baseline fit: %7.2f km/sec ' % (minvel))
        print('20 channels near this velocity will be used for the fit')
    elif sys.argv[iarg].upper() == '-VB':   # now look for flags with arguments
        iarg = iarg+1
        maxvel = float(sys.argv[iarg])
        print('Maximum velocity for baseline fit: %7.2f km/sec ' % (maxvel))
        print('20 channels near this velocity will be used for the fit')
    else:
        break
    iarg = iarg + 1
    timearg = iarg
    namearg = timearg+1
# end of while not reading file names

# first argument is average time
try:
    avetimesec = float(sys.argv[timearg])
except:
    print('First position argumen must be the average time in seconds')
    print('First argument: ', sys.argv[timearg])
    exit()

avetime = dt.timedelta(seconds=avetimesec)  # need datatime format

if doDebug:
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

tmin = 20.0  # Tsys never less than 20 K
tmax = 999.0 # define reasoanable value limits

velbaseline = [ minvel, maxvel] # array of velocities in km/sec

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
firstaz = -1
otheraz = -1
minGlat = +90.
maxGlat = -90.
minel = 200.
maxel = -200.

DT = dt.timedelta(seconds=0.)

#rest of arguments are the file names
names = sys.argv[namearg:]
names = sorted(names)
nFiles = len(names)

def fit_baseline( xs, ys, imin, imax, nchan):
    """
    fit baseline does a polynomical fit over channels in a select range
    The baseline is returned.
    """
    
    xfit = np.concatenate( (xs[imin-nchan:imin+nchan],xs[imax-nchan:imax+nchan]))
    yfit = np.concatenate( (ys[imin-nchan:imin+nchan],ys[imax-nchan:imax+nchan]))
# calculate polynomial
    z = np.polyfit(xfit, yfit, 3)
    f = np.poly1d(z)

# calculate new x's and y's
    yout = f(xs)

    return yout

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
        tsys[icenter] = (tsys[icenter-2] + tsys[icenter+2])*.5
        tsys[icenter-1] = (3.*tsys[icenter-2] + tsys[icenter+2])*.25
        tsys[icenter+1] = (tsys[icenter-2] + 3.*tsys[icenter+2])*.25

    # get indicies for given max and min velocities
    imin, imax = velocity_to_indicies( vel, minvel, maxvel)

    if doPoly:
        baseline = fit_baseline( vel, tsys, imin, imax, nfit)
        tsys = tsys - baseline
    else:
        ya = np.median(tsys[(imin-nFit):(imin+nFit)])
        yb = np.median(tsys[(imax-nFit):(imax+nFit)])
        slope = (yb-ya)/np.float(imax-imin)
#                baseline = tsys
#                print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
# This is the only difference between M and T processing
        if doDebug:
            print("Slope %7.3f K/channel" % (slope))
        for iii in range( nData):
            tsys[iii] = tsys[iii] - (ya + (slope*(iii-imin)))
        
    return tsys

if doDebug:
    print('First Name: ',names[0], namearg)

hotnames = copy.deepcopy( names)
coldnames = copy.deepcopy( names)
plotnames = copy.deepcopy( names)

############################################################ Start Hot processing
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

#previously just copy
hv = copy.deepcopy(yv)
if flagRfi:
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

vel = np.zeros(nData)
# create index array
for jjj in range (nData):
    vel[jjj] = c * (nureference - xv[jjj])/nureference

xa, xb = velocity_to_indicies( vel, minvel, maxvel)
print('Min Vel, channel: ',xa, minvel)
print('Max Vel, channel: ',xb, maxvel)
                                   
if doDebug:
    print('Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % (minGlat, maxGlat))

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

# if hanning smoothing the hot data
if doHanning:
    yv = copy.deepcopy(hv)
    for jjj in range(nData-2):
        iii = jjj+1
        hv[iii] = (yv[iii-1]+yv[iii]+yv[iii]+yv[iii+1])*.25

if doDebug:
    print("hv: %8.2f (counts)" % (hv[300]))

if doSave:
    hot.writecount = 1
    hot.ydataA = hv
    hot.write_ascii_ave( outDir)
else:
    print("")
        
if doDebug: 
    print('')
    print('Hot File Coord: ')
    print('nChan, refChan: ',hot.nChan, hot.refChan)
    print('centerFreqHz  : ',hot.centerFreqHz)
    print('bandwidthHz   : ',hot.bandwidthHz)
    print('')
    print('N Channels    : ',nData)
    print('velocity range: ', velbaseline)
    print('Reference Freq. ', nureference)
    print('Min Vel at channel: ',xa, minvel)
    print('Max Vel at channel: ',xb, maxvel)
                                   

############################################################ Start Cold processing
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

cv = copy.deepcopy(yv)
if flagRfi:
    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

# finally compute gain on a channel by channel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)
    vel[iii] = c * (nureference - xv[iii])/nureference

#now want gain using only hot counts, but must add in Tsys
tSys = cv/gain
tSysMiddle = np.median(tSys[n6:n56])
CountsPerKelvin = np.median(gain[n6:n56])

# for remainder of calculations only use hot counts for calibration
for iii in range(nData):
    gain[iii] = hv[iii]/(thot + tSysMiddle - tcold)
    
if doDebug:
    print("cv: %8.2f (counts)" % (cv[300]))

if doSave:
    high.writecount = 2
    high.ydataA = copy.deepcopy(cv)
    high.write_ascii_ave( outDir)
else:
    print("")

########################### End of Cold Spectrum
fig, ax1 = plt.subplots(figsize=(10, 6))

if not doCalibrate: # if not calibrating, show hot and cold spectra
    tsky = copy.deepcopy(hv)
    if doSub:
        ya = np.median(tsky[(xa-nFit):(xa+nFit)])
        yb = np.median(tsky[(xb-nFit):(xb+nFit)])
        slope = (yb-ya)/(xb-xa)
        for iii in range( nData):
            tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))
    yallmin = min(tsky[xa:xb])
    yallmax = max(tsky[xa:xb])
    plt.plot(vel, tsky, '-r', linestyle='-', label="Hot  Load", lw=2)
    if doDebug:
        print('hv[300]: ', tsky[300])

    tsky = copy.deepcopy(cv)
    if doSub:
        ya = np.median(tsky[(xa-nFit):(xa+nFit)])
        yb = np.median(tsky[(xb-nFit):(xb+nFit)])
        slope = (yb-ya)/(xb-xa)
        for iii in range( nData):
            tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))
    # reset the min/max values to include the hot and cold loads
    ymin = min(tsky[xa:xb])
    ymax = max(tsky[xa:xb])
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)

    plt.plot(vel, tsky, '-b', linestyle='-', label="Cold Load", lw=2)
    if doDebug:
        print('cv[300]: ', tsky[300])

# end if not calibrating, then plot hot/cold

nRead = 0        
# create structure to hold spectrum
rs = radioastronomy.Spectrum()

# now read through all data and average cold sky obs
for filename in plotnames:
    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    extension = parts[1]
    nRead = nRead + 1
# exclude hot load data for averaging
    if nRead == nFiles:   # if this is the last file, must force output
        allFiles = True

# also exclude summaries
    if extension == 'sum':
        continue

#  print filename
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if rs.telel < 0. and (not allFiles):
        continue

    if doFold:
        rs.foldfrequency()

    date, time = rs.datetime()

    if firstdate == "":            # if very first spectrum
        firstdate = date
        firstaz = rs.telaz
        lastaz = rs.telaz          # heep track of a range of azimuths
        otheraz = rs.telaz
        veryfirstutc = rs.utc
        verylastutc = rs.utc

    if otheraz != rs.telaz:
        otheraz = rs.telaz

    if minel > rs.telel:
        minel = rs.telel
    if maxel < rs.telel:
        maxel = rs.telel

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
        firstutc = rs.utc
        lastutc = rs.utc
        ncold = 0

    newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
    newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl
    if newObs:
        if lastfreq != rs.centerFreqHz:
            print("Change LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz")
            lastfreq = rs.centerFreqHz
        if lastbw != rs.bandwidthHz:
            print("Change LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz")
            lastbw = rs.bandwidthHz
        if lastgain != rs.gains[0]:
            print("Change LastGain: ", lastgain, "New: ", rs.gains[0], " dB")
            lastgain = rs.gains[0]
        if newAzEl:
            print("Change LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees")
            lastaz = rs.telaz
            lastel = rs.telel

    # must add this file in if it is the last file and not a new file
    if allFiles and (not newObs):
        cold, rs, ncold, firstutc, lastutc = average_spec( cold, rs, ncold, firstutc, lastutc)

    DT = dt.timedelta(seconds=rs.durationSec)

    if ncold > 0 and (not newObs):
        # time difference is between mid-points of integrations
        DT = DT + (rs.utc - cold.utc)
        lastdate = date
        if rs.utc > verylastutc:
            verylastutc = rs.utc
        if rs.utc < veryfirstutc:
            veryfirstutc = rs.utc

    if doDebug:
        print("Ncold: ",ncold)
        print("newObs: ",newObs)
        print("newObs: ",allFiles)
        
    # if time to average (or end of all files)
    if (DT > avetime) or newObs or allFiles:
        cold.ydataA = cold.ydataA/float(cold.durationSec)
        if doDebug:
            print("Average duration: %7.1f " % (cold.durationSec))

        yv = cold.ydataA
        if flagRfi:
            yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

        # compute velocity correction for this direction and date
        corr = compute_vbarycenter( cold)
        velcorr = vel + corr
        xa, xb = velocity_to_indicies( velcorr, minvel, maxvel)

        xmin = min(xv[xa:xb])
        xmax = max(xv[xa:xb])
        xallmin = min(xmin, xallmin)
        xallmax = max(xmax, xallmax)
        count = cold.count
        note = cold.noteA

        ncolor = min(nmax-1, nplot)  # select color for this spectrum
        tsky = compute_tsky_hotcold( yv, gain, velcorr, minvel, maxvel)
            
        ymin = min(tsky[xa:xb])
        ymax = max(tsky[xa:xb])
        yallmin = min(ymin,yallmin)
        yallmax = max(ymax,yallmax)
        # compute average time from first and last utcs
        aveutc,duration = radioastronomy.aveutcs( firstutc, lastutc)
        if doDebug and (firstutc == lastutc):
            print("first and last utcs are the same: ", firstutc)
        # keep the average time, but use the duration from the integration times.
        cold.utc = aveutc 
            # pull out coordinates for labeling
        az = cold.telaz
        el = cold.telel
        cold.azel2radec()    # compute ra,dec from az,el and average utc
        gallon = cold.gallon
        gallat = cold.gallat

        tSys = np.median(tsky[xa:xb])
        tStdA = np.std(tsky[(xa-nFit):xa])
        tStdB = np.std(tsky[xb:(xb+nFit)])
        tStd = (tStdA+tStdB)/.2

        if doDebug:
            print(('Calibrated Max, Median and Std Dev: %8.3f %8.3f %8.3f (K)' % (ymax, tSys, tStd)))
            print('xa,xb: ',xa,xb)
        date, time = cold.datetime()
        labeltime = date + " " + time
        if minel == maxel:
            label = '%s L,L=%5.1f,%4.1f' % (labeltime, cold.gallon, cold.gallat)
        else:
            label = '%s A,E=%.0f,%4.0f L,L=%5.1f,%4.1f' % (labeltime, cold.telaz, cold.telel, cold.gallon, cold.gallat)
        outstring = ' Max: %9.1f  Median: %8.2f +/- %5.2f %s %3d %s' % (ymax, tSys, tStd, cold.bunit[0], ncold, label)
        sys.stdout.write(outstring)
        sys.stdout.flush()

        lw = 1
        if cold.gallat < 7.5 and cold.gallat > -7.5:
            lw=4
        elif cold.gallat < 15. and cold.gallat > -15.:
            lw=2
        else:
            lw=1
        plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=lw)

        nplot = nplot + 1
        if doKelvins:
            cold.writecount = nplot+2
            cold.ydataA = tsky
            cold.write_ascii_ave( outDir)
        else:
            print("")
                
        # indicate that we're now restarting averaging
        ncold = 0
        # end if a new observation

    if allFiles:
        break
    
    # Average in most recently read spectrum
    cold, rs, ncold, firstutc, lastutc = average_spec( cold, rs, ncold, firstutc, lastutc)

    if ncold == 1:
        # At each sum restart, must recompute velocity axis
        xv = cold.xdata
        vel = copy.deepcopy( xv)
        vel = nureference - vel         # converting to frequency offset
        vel = (c/nureference) * vel     # now conver to velocity
        lastutc = cold.utc
    # end if not a enough time

    #end for all files to sum
# end of all files to read

# if observations cover multiple days
if firstdate != lastdate:
    date = firstdate + " - " + lastdate
else:
    date = firstdate

if mytitle == "":
    mytitle = rs.site

if (firstaz == otheraz):
    mytitle = mytitle + "; Az = %6.1f, " % (firstaz)
else:
    mytitle = mytitle + "; Az = %6.1f to %6.1f, " % (firstaz, otheraz)

if minel == maxel:
    mytitle = mytitle + " El=%6.1f" % (minel)
else:
    mytitle = mytitle + " El=%6.1f to %6.1f" % (minel, maxel)

fig.canvas.set_window_title(mytitle)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_rotation('vertical')
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_rotation('vertical')
# select the relevant range of velocities (km/sec) for plotting
plt.xlim(minvel, maxvel)
# keep the plot from becoming too narrow
if yallmax < 8:
    yallmax = 8
# set the y scale to go above and below all data
plt.ylim((yallmin*.97)-1., 1.15*yallmax)
plt.title(mytitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
if doCalibrate:
    plt.ylabel('Intensity (Kelvins)', fontsize=16)
else:
    plt.ylabel('Intensity (Counts)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
