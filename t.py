#Python Script to plot calibrated  NSF spectral integration data.
#plot the raw data from the observation
#HISTORY
#20APR25 GIL use average_spec() and normalize_spec() for hot+cold
#20APR24 GIL optionally save raw hot and cold load obs
#20APR22 GIL improve estimate of integrated intensity
#20APR21 GIL fix up integrated intensities and velocities
#20APR17 GIL add option to plot baseline-subtracted obs
#20APR16 GIL complete logging of calibration values
#20APR14 GIL add log of calibration values
#20APR01 GIL reduce printout
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
import gainfactor as gf

try:
    from PyAstronomy import pyasl
    baryCenterAvailable = True
except:
    print("!!!! Unable to import PyAstronomy !!!!")
    print("Can not compute Bary Center velocity offset")
    print("In Linux, try: ")
    print("sudo pip install PyAstronomy")
    print("or")
    print("sudo pip3 install PyAstronomy")
    baryCenterAvailable = False

# Maximum number of Plots
maxPlot = int(25)
# default values
avetimesec = 3600.
# put your list of known RFI features here.  Must have at least two.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]
# min and maximum default velocities
maxvel = 220.
minvel = -maxvel
# min and max velocities for intensity integration
maxSVel = 150.  # symetric integration
minSVel = -maxSVel
# cpu Index for normalizing gain values
cpuIndex = 0
# Keep/plot baseline Subtracted spectra
doBaseline = False

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
# writeTsys files
writeTsys = True
writeTsys = False
# to address an old problem, optionally allow folding spectra
doFold = False
# specify lowest  elevation for cold load averge
lowel = 60.
# define fitOrder for intensity measurement
fitOrder = int(1)

# optionally turn on debug plotting
doDebug = False
myTitle = ""      # Default no input title
saveFile = ""     # Default no saveFileName
hotLoadFile = ""
coldLoadFile = ""

iarg = 1
if nargs < 3:
    print("T: Comput Tsys calibrated horn observations")
    print("Usage: T [-F <order>] [-L <velocity>] [-H <velocity>] <average_seconds> <files>")
    print("Where <average_seconds>: Number of seconds of observations to average.")
    print("-B optionally plot/keep the baseline subtratcted spectra")
    print("-C optionally flag the center of the band")
    print("-F optionally do a polynomial baseline fit")
    print("-L optionally set the low velocity region for baseline fit")
    print("-I optionall set Processor Index")
    print("-H optionall set the high velocity region for baseline fit")
    print("-R optionally flag known RFI lines")
    print("-S <filename> optionally set summary file name")
    print("-W optionally write the calibrated Tsys files")
    print("-MINEL optionally set the lowest elevation allowed for calibration obs (default 60d)")
    print("Observation file list must include at least one hot load file")
    print("")
    print("Glen Langston - NSF   November 22, 2019")
    exit()

# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-B':
        doBaseline = True
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-D':
        print( 'Adding Debug Printing')
        doDebug = True
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
    elif sys.argv[iarg].upper() == '-G':
        iarg = iarg+1
        hotLoadFile = sys.argv[iarg]
        iarg = iarg+1
        coldLoadFile = sys.argv[iarg]
        print("Gain calibration Using Hot Load Obs: %s and Cold Obs: %s" % (hotLoadFile, coldLoadFile))
    elif sys.argv[iarg].upper() == '-H':
        iarg = iarg+1
        maxvel = np.float( sys.argv[iarg])
        print('Maximum (high) velocity for sum: %7.2f km/sec' % (maxvel))
    elif sys.argv[iarg].upper() == '-I':
        iarg = iarg+1
        cpuIndex = np.float( sys.argv[iarg])
        print('Telescope Index: %2d' % (cpuIndex))
    elif sys.argv[iarg].upper() == '-L':
        iarg = iarg+1
        minvel = np.float( sys.argv[iarg])
        print('Minium (low)  velocity for sum: %7.2f km/sec' % (minvel))
    elif sys.argv[iarg].upper() == '-N':   # if number of spectra to plot
        iarg = iarg+1
        maxPlot = int(sys.argv[iarg])
        if maxPlot < 1:
            print("Not Plotting")
        else:
            print("Plot will have a maximum of %d spectra" % (maxPlot))
    elif sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-T':   # if plot title provided
        iarg = iarg+1
        myTitle = sys.argv[iarg]
        print( 'Plot Title : ', myTitle)
    elif sys.argv[iarg].upper() == '-S':   # if save file name provided
        iarg = iarg+1
        saveFile = sys.argv[iarg]
        print( 'Save File: ', saveFile)
    elif sys.argv[iarg].upper() == '-MINEL':  # if min elevation 
        iarg = iarg + 1
        lowel = float( sys.argv[iarg])
        print( "Using elevations > %7.2f (d) for Cold load calculation" % (lowel))
    elif sys.argv[iarg].upper() == '-W':
        writeTsys = True
    elif sys.argv[iarg].upper() == '-Z':
        doFold = True
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
xallmax = -9.e9
xallmin = 9.e9
ymin = 9.e9
ymax = -9.e9
yallmax = ymax
yallmin = ymin
# velocities for fitting baselines

c = 299792.458  # (Speed of light  km/sec)
nuh1 = 1420.40575 # neutral hydrogen frequency (MHz)
thot = 285.0  # define hot and cold load temperatures
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
ave_hot  = radioastronomy.Spectrum()
ave_cold = radioastronomy.Spectrum()

def average_spec( ave_spec, in_spec, nave, firstutc, lastutc):
    """
    Averages spectra and deal with issues of weighting functions by observing duration
    input/output
    ave_spec   spectrum structure containing weighted average
    input:
    in_spec    spectrum to be added to the average
    in/out:
    nave       Count of number of spectra averaged so far.
               nave = 0 implies initializing the sum.
    firstutc   date of first observation averaged
    lastutc    date of current last date of spectrum to be averaged
    """

    if doDebug:
        medianData = np.median( in_spec.ydataA[n6:n56])
    # remove number of spectra averaged scaling.
    in_spec.ydataA = (in_spec.ydataA/in_spec.count) 
    if doDebug:
        medianScale = np.median( in_spec.ydataA[n6:n56])
        print( "Input: %8.3f, count: %6d; Scale: %8.3f" % (medianData, in_spec.count, medianScale))

    # if restarting the sum
    if nave == 0:
        ave_spec = copy.deepcopy(in_spec)  # initial spectrum is one just read
        firstutc = in_spec.utc
        lastutc = in_spec.utc
        nave = 1
        # print( 'Xmin: ', min(ave_spec.xdata)/1e6, 'Xmax: ', max(ave_spec.xdata),' MHz')
        # sums are weighted by durations
        ave_spec.ydataA = in_spec.ydataA * in_spec.durationSec  # replace with duration scaleing
        # keep track of observing time for weighted sum
        ave_spec.durationSec = in_spec.durationSec
    else: # else not enough time yet, average ave_spec data
        if in_spec.utc < firstutc:
            firstutc = in_spec.utc
        elif in_spec.utc > lastutc:
            lastutc = in_spec.utc
        ave_spec.count = ave_spec.count + in_spec.count
        nave = nave + 1
        ave_spec.ydataA = ave_spec.ydataA + (in_spec.ydataA * in_spec.durationSec)
        # keep track of observing time for weighted sum
        ave_spec.durationSec = ave_spec.durationSec + in_spec.durationSec
    return ave_spec, in_spec, nave, firstutc, lastutc

def normalize_spec( ave_spec, nave, firstutc, lastutc):
    """
    Normaize the average after completing sum
    input/output
    ave_spec   input raw sum of observation 
               output normalized sum of observations
    nave       input number of spectra averaged, output zero
    input      first and last utcs in observation
    Note:  The count of observations is not used, rather the integration
    is weighted by durations.  This corrects for observations with 
    different integration times.
    """
    # now renormalize for total integration time
    ave_spec.ydataA = ave_spec.ydataA/float(ave_spec.durationSec)
    # compute average time from first and last utcs
    aveutc, duration = radioastronomy.aveutcs( firstutc, lastutc)
    ave_spec.utc = aveutc
    #need to re-calculate representative RA,Dec for average time
    ave_spec.azel2radec()
    nave = 0
    return ave_spec, nave

# record the time range of hot load observations and count of observations
firstutc = 0
lastutc = 0    
nhot = 0

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

    # if a hot load observation
    if rs.telel < 0:
        if doFold:
            rs.foldfrequency()
            
        # accumulate spectra
        ave_hot, rs, nhot, firstutc, lastutc = \
            average_spec( ave_hot, rs, nhot, firstutc, lastutc)
    else: # else above horizon, find min, max galactic latitudes
        if rs.gallat > maxGlat:
            maxGlat = rs.gallat
        if rs.gallat < minGlat:
            minGlat = rs.gallat
        if minel > rs.telel:
            minel = rs.telel
        if maxel < rs.telel:
            maxel = rs.telel
        # track if multiple azimuths
        if firstaz < 0:
            firstaz = rs.telaz
            otheraz = firstaz
        if firstaz != rs.telaz:
            otheraz = rs.telaz

    # end for all files

if nhot > 0:
    print( "Found %3d Hot load observations" % (nhot))
    ave_hot, nhot = normalize_spec( ave_hot, nhot, firstutc, lastutc)
else:
    print( "No Hot load data, can not calibrate")
    exit()

# convert to MHz
xv = ave_hot.xdata * 1.E-6
hv = copy.deepcopy(ave_hot.ydataA)

#previously just copy
if flagRfi:
    yv = copy.deepcopy(hv)
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
    ave_hot.ydataA = hv

nData = len( xv)        # get data length and indicies for middle of spectra 
n6 = int(nData/6)
n26 = int(2*n6)
n46 = int(4*n6)
n56 = int(5*n6)

vel = np.zeros(nData)
# create index array
for jjj in range (0, nData):
    vel[jjj] = c * (nuh1 - xv[jjj])/nuh1

xa, xb = gf.velocity_to_indicies( vel, minvel, maxvel)

if doDebug:
    print( 'Min Vel  %7.1f, Max Vel  %7.1f' % ( minvel, maxvel))
    print( 'Min Chan %7d, Max Chan %7d' % (xa, xb))
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
#    print( "Spectrum count: %d " % (rs.count))
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

    if rs.telel < lowel:  #if elevation too low for a cold load obs
        continue

    if rs.gallat > maxGlat or rs.gallat < minGlat:
        if nhigh == 0:
            cold = copy.deepcopy( rs)
            cold.ydataA = (rs.ydataA/rs.count) * rs.durationSec
            cold.count = rs.count
            cold.durationSec = rs.durationSec
            nhigh = 1
        else:
            cold.ydataA = cold.ydataA + ((rs.ydataA/rs.count) * rs.durationSec)
            cold.count = cold.count + rs.count
            cold.durationSec = cold.durationSec + rs.durationSec
            nhigh = nhigh + 1

if nhigh < 1.:
    print( "No high galactic latitude data: can not calibrate")
    exit()
else:
    cold.ydataA = cold.ydataA/cold.durationSec
    print( "Found %3d High Galactic Latitude spectra" % (nhigh))

yv = copy.deepcopy(cold.ydataA)
cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

# compute gain on a channel by channel basis for tRx calculation
gainHC = np.zeros(nData)
for iii in range(nData):
    gainHC[iii] = (hv[iii] - cv[iii])/(thot - tcold)

trx = np.zeros(nData)
for iii in range(nData):
    trx[iii] = (cv[iii]/gainHC[iii]) - tcold

#now prepare to compute tRx, which is based only on cold load observations
tRxA = np.median(trx[n6:n26])
tRxB = np.median(trx[n46:n56])
tRxMiddle = (tRxA + tRxB)*.5

tStdA = np.std(trx[n6:n26])
tStdB = np.std(trx[n46:n56])
tRms  = (tStdA + tStdB) * .5

print( "Median Receiver Temp: %7.2f +/- %5.2f (%5.2f %5.2f) (K)" % ( tRxMiddle, tRms, tStdA, tStdB))

# for remainder of calculations only use hot counts for calibration
# Using hot load only reduces interference effects
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = hv[iii]/(thot + tRxMiddle)

gainA = np.median(gain[n6:n26])
gainB = np.median(gain[n46:n56])
gainAve = 2.0/(gainA + gainB)  # Report gain in K per Count 

def compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold):
    """
    Compute TSky based on hot and coold load observations
    Note that here the yv[], hv[] and cv[] values are normalized by the total number of integrations.
    The raw data, read in, are increased linearly by the number of spectra averaged.
    """
    nData = len(xv)
    vel = np.zeros(nData)
    tsky  = np.zeros(nData)    # initialize arrays

    for jjj in range (0, nData):
        vel[jjj] = c * (nuh1 - (xv[jjj]))/nuh1
        tsky[jjj] = yv[jjj]/gain[jjj]

    if flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
        icenter = int(nData/2)
        tsky[icenter] = (tsky[icenter-2] + tsky[icenter+2])*.5
        tsky[icenter-1] = (3.*tsky[icenter-2] + tsky[icenter+2])*.25
        tsky[icenter+1] = (tsky[icenter-2] + 3.*tsky[icenter+2])*.25

    return tsky, vel
# end of compute_tsky_hotcold()

# if plotting
if maxPlot > 0:
    fig, ax1 = plt.subplots(figsize=(10, 6))

# compute the reciever temperature
trx = np.zeros(nData)
zeros = np.zeros(nData)

avetime = datetime.timedelta(seconds=avetimesec)

nRead = 0        
nave = 0

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

#   create spectrum structure, empty
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
        ave_spec = copy.deepcopy( rs)
        firstutc = rs.utc
        lastutc = rs.utc
        nave = 0

    if nave > 0:
        # time difference is between mid-points of integrations
        dt = rs.utc - ave_spec.utc 
        # add the time since midpoint of latests
        dt = dt + datetime.timedelta(seconds=rs.durationSec)
        lastdate = date

        newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
        newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl
        # if this is the last file and there was not a change in observing setup 
        if allFiles and (not newObs):
            if rs.nChan != ave_spec.nChan:
                print("File Size Error: %d != %d for file %s; Skipping ..." % (rs.nChan, ave_spec.nChan, filename))
            else:
                cold, rs, nave, firstutc, lastutc = average_spec( cold, rs, nave, firstutc, lastutc)

        # if time to average (or end of all files)
        if (dt > avetime) or newObs or allFiles:
            if doDebug:
                medianData = np.median( ave_spec.ydataA[n6:n56])
                print( "Average duration: %7.1f, Median:  %8.3f" % (ave_spec.durationSec, medianData))

            # not calibrating hot load observations.
            if ave_spec.telel < 0.:
                # Reset the ncold count and restart sum
                nave = 0
                continue
            ave_spec, nave = normalize_spec( ave_spec, nave, firstutc, lastutc)
            xv = ave_spec.xdata * 1.E-6  # covert to MHz
            yv = ave_spec.ydataA 
            if flagRfi:
                yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

            xmin = min(xv)
            xmax = max(xv)
            xallmin = min(xmin, xallmin)
            xallmax = max(xmax, xallmax)
            count = ave_spec.count
            note = ave_spec.noteA
                    #print('%s' % note)
            ncolor = min(nmax-1, nplot) 

            tsky, vel = compute_tsky_hotcold( xv, yv, hv, cv, thot, tcold)
            # get tsys from averages of ends of spectra
            tSys = np.median(tsky[n6:n56])
            tStdA = np.std(tsky[n6:n26])
            tStdB = np.std(tsky[n46:n56])
            cA = np.median(cv[n6:n26])
            cB = np.median(cv[n46:n56])
            counts = (cA+cB)/2.
            tStd = (tStdA+tStdB)/2.
            ave_spec.tSys = tSys
            ave_spec.tRx = tRxMiddle
            ave_spec.tRms = tStd
            ave_spec.tint = ave_spec.durationSec
            ave_spec.bunit = 'Kelvins'
            ave_spec.KperC = gainAve
            
            # compute velocity correction for this direction and date
            corr = gf.compute_vbarycenter( cold)
            velcorr = vel + corr

            # compute indicies for min and max velocity
            imin, imax = gf.velocity_to_indicies( velcorr, minvel, maxvel)

            # compute indicies for min and max velocity to integrate
            iVmin, iVmax = gf.velocity_to_indicies( velcorr, minSVel, maxSVel)


            # keep the average time, but use the duration from the integration times.
            # pull out coordinates for labeling
            az = ave_spec.telaz
            el = ave_spec.telel
            ave_spec.azel2radec()    # compute ra,dec from az,el and average utc
            gallon = ave_spec.gallon
            gallat = ave_spec.gallat
            label = 'L,L=%5.1f,%6.1f' % (gallon, gallat)
            
            # this code computes and subtracts a baseline so that source intensities can be compared.
            # next compute the integrated intensities after baseline subtraction
            baseline = gf.fit_baseline( velcorr[0:nData], tsky[0:nData], imin, imax, 10, fitOrder)

            # remove baseline to get the source spectrum
            tSource = tsky[0:nData] - baseline[0:nData]

            # if plotting/keeping the baseline subtracted spectra, transfer to Sky
            if doBaseline:
                tsky = tSource

            # set min and max y for plotting (only)
            ymin = min(tsky[imin:imax])
            ymax = max(tsky[imin:imax])

            # now compute integrated intensity and noise estimates
            nv = iVmax - iVmin
            # create sub-arrays of intensity and velocity
            tSs = tSource[iVmin:iVmax]
            vSs = velcorr[iVmin:iVmax]
            tSourcemin = min(tSs)
            tSourcemax = max(tSs)
            # get index to maximum value; then get velocity
            iSourcemax = np.argmax(tSs)
            velSource = vSs[iSourcemax]

            # integrate over spectrum for required velocity range
            tSum = np.sum(tSs)
#            tSumRms = np.std(tSs)*np.sqrt(nv)
            tSumRms = np.std(tSs)
            # Integration is reported in Kelvin*Km/Sec;  Multiply by source velocity range
            tSumKmSec = tSum * ( maxSVel - minSVel)
            dTSumKmSec = tSumRms * ( maxSVel - minSVel)

            # computed the integrated velocity
            tvs = tSs*vSs
            tVSum = np.sum(tvs)
            tVSumRms = np.std(tvs)

            # now compute intensity weighted velocity
            if tSum > 0.:
                tVSum = tVSum/tSum
                tVSumRms = tVSumRms/tSum
            else:
                tVSum=0.

            nChan2 = int(ave_spec.nChan/2)
            dV = (velcorr[nChan2] - velcorr[nChan2 - 4])*.25
            if dV < 0:
                dV = - dV
            # print diagnostic of integration
            if (doDebug):
                print("Min, max Velocity   : %7.1f  to %6.1f; %d,%d" % (minvel,maxvel, imin, imax))
                print("Min, max Velocity I : %7.1f  to %6.1f; %d,%d" % (minSVel,maxSVel, iVmin, iVmax))
                print("Average Intensity(K): %7.3f +/- %6.3f K (%d)" % (tSum/nv, tSumRms/nv, nv))
                print("Int.  Vel.  K km/s  : %7.0f +/- %6.0f K*km/sec" % (tSumKmSec,  dTSumKmSec))
                print("Peak    Velocity    : %7.1f +/- %6.1f km/sec (%d)" % (velSource,  dV, iSourcemax))
                print("Integrated Velocity : %7.1f +/- %6.1f km/sec" % (tVSum,  tVSumRms))

            avedatetime = ave_spec.utc.isoformat()
            datestr = avedatetime.split('T')
            atime = datestr[1]
            timeparts = atime.split('.')
            labeltime = timeparts[0]
            label = '%s, A,E: %5s,%5s, L,L: %5.1f,%6.1f' % (labeltime, az, el, gallon, gallat)
            if minel == maxel:
                label = '%s L,L=%5.1f,%5.1f (%d)' % (labeltime, gallon, gallat, nave)
            else:
                label = '%s L,L=%5.1f,%5.1f A,E=%4.0f,%4.0f' % (labeltime, gallon, gallat, az, el)
            print( ' Max: %9.1f  Median: %8.2f +/- %5.2f %3d %s' % (tSourcemax, tSys, tStd, nave, label))
            if int(nplot) < int(maxPlot):
                yallmin = min(ymin,yallmin)
                yallmax = max(ymax,yallmax)
                if gallat < 7.5 and gallat > -7.5:
                    plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
                elif gallat < 15. and gallat > -15.:
                    plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
                else:
                    plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            nplot = nplot + 1
            # this code computes and subtracts a baseline so that source intensities can be compared.

            gf.saveTsysValues( saveFile, cold, cpuIndex, tSourcemax, velSource, dV, tVSum, tVSumRms, tSumKmSec, dTSumKmSec)
            if writeTsys:
                ave_spec.ydataA = tsky
                outname = radioastronomy.utcToName( aveutc)
                outname = outname + ".kel"  # output in Kelvins
                ave_spec.count = 1
                ave_spec.write_ascii_file("../", outname)                
# flag restarting the sum
            nave = 0

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


    # end if a cold file
    if allFiles:
        break

    if rs.telel > 0:
    # Average in most recently read spectrum
        if rs.nChan != ave_spec.nChan:
            print("File Size Error: %d != %d for file %s; Skipping ..." % (rs.nChan, ave_spec.nChan, filename))
        else:
            ave_spec, rs, nave, firstutc, lastutc = average_spec( ave_spec, rs, nave, firstutc, lastutc)

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
    myTitle = "%s   " % (date)
else:
    myTitle = myTitle + "  "


if cpuIndex > 0:
    myTitle = myTitle + ("T:%d " % (cpuIndex))

if (firstaz == otheraz):
    myTitle = myTitle + "Az = %6.1f, " % (firstaz)
else:
    myTitle = myTitle + "Az = %6.1f to %6.1f, " % (firstaz, otheraz)

if minel == maxel:
    myTitle = myTitle + " El=%6.1f" % (minel)
else:
    myTitle = myTitle + " El=%6.1f to %6.1f" % (minel, maxel)

if (maxPlot < 1) or (nplot < 1):
    exit()
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
#plt.xlim(-250., 250.)
plt.xlim(minvel, maxvel)
# keep the plot from becoming too narrow
dy = yallmax - yallmin
if dy < 8:
    dy = 8

if doDebug:
    print ( "Y min, max: %8.3f, %8.3f" % (yallmin, yallmax))

# set the y scale to go above and below all data
plt.ylim((yallmin-(dy/8.)), (yallmax+(dy/4.)))
plt.title(myTitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.show()
