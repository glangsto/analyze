"""
t.py computes the calibrated temperature versus frequency or velocity
spectrum of radio astronomy observations.
Glen Langston National Scioence Foundation
"""
#Python Script to plot calibrated  NSF spectral integration data.
#plot the raw data from the observation
#HISTORY
#25APR28 GIL update for labels without elevation change
#24OCT22 GIL fix font size for labels
#23SEP27 GIL recalculate channel ranges for baseline subtraction for each ave.
#23Jul05 GIL use galactic hydrogen observations above 40 degrees
#23Apr25 GIL use cold load for velocity computations
#23Apr21 GIL debug spectral line obs at 1612 MHz
#22Jun11 GIL Write .kel (vins) file with channel/velocity/freq units of plot
#22Apr28 GIL fix use of cv array before initialization
#22Apr28 GIL debug median filtering
#22APR22 GIL put some functions in hotcold.py
#22APR20 GIL find channels at the velocities identified
#22MAR27 GIL median option for output vector
#22FEB18 GIL Bucket Horn RFI lines
#22FEB02 GIL add additional RFI lines temporarily
#21JUN26 GIL add option to write velocities
#21MAY09 GIL optionally flag 1422.0 MHz
#21FEB23 GIL add a zero intensity line option
#21JAN19 GIL fix normalizations again
#21JAN15 GIL add a little more documentation
#21JAN06 GIL fix normalizations
#20DEC28 GIL check for a average time in first argument
#20DEC16 GIL flag and ignore short spectra files
#20NOV12 GIL fix plotting to files
#20AUG28 GIL add option to only write the plots
#20JUL27 GIL add options to set thot, tcold, fix writing Kelvins
#20JUN02 GIL plot frequency if velocity out of range
#20MAY06 GIL update help
#20APR30 GIL update help
#20APR28 GIL write hot and cold obs
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
#18FEB06 GIL keep track of first az,el
#17AUG18 GIL plot the last spectrum in the series
#17AUG17 GIL Note elevation range
#16Oct07 GIL check for changes in frequency and/or gain
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import sys
import os
import copy
import datetime
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import gainfactor as gf
import radioastronomy
import interpolate
import tsys
import hotcold

# Maximum number of Plots
maxPlot = int(25)
fileTag = ""

# default values
avetimesec = 3600.
# put your list of known RFI features here.  Must have at least two.
#linelist = [1400.00, 1418.21, 1420.0, 1422.0]  # RFI lines in MHz
linelist = [1400.00, 1419.92, 1420.0, 1421.24, 1421.39, 1421.54]  # RFI lines in MHz
linelist = [1400.00, 1419.49, 1419.63, 1419.93, 1420.0, 1420.37, 1420.66, 1420.95, 1421.1, 1421.25, 1421.39, 1421.54]  # RFI lines in MHz
linewidth = [5, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
nlist = len (linelist)
nwidth = len(linewidth)
# min and maximum default velocities
maxvel = 175.
minvel = -maxvel
# min and max velocities for intensity integration
maxSVel = 150.  # symetric integration
minSVel = -maxSVel
# cpu Index for normalizing gain values
cpuIndex = 0
# Keep/plot baseline Subtracted spectra
doBaseline = False

#first argument is the averaging time in seconds
timearg = 1
namearg = 2

# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False
# put list of RFI features here, for interpolation later
flagRfi = False
# write log of integrated intensities
writeTsys = False
# write calibrated intensities, Kelvin
writeKelvin = False
# Usually writing file header, but optionally skip
doWriteHeader = True
# if limiting velocity output to specified min and max velocity
doLimitVel = False
# if plotting to a file, specify the directory
doPlotFile = False
plotFileDir = "~/"
# if keeping average hot and cold load
doKeep = False
# specify lowest elevation for cold load averge
lowel = 10.
# define fitOrder for intensity measurement
fitOrder = int(1)
# define flag to compute barycentric velocity
doBary = True
minBaryVel = 0.
maxBaryVel = 0.
# optionally turn on debug plotting
doDebug = False
myTitle = ""      # Default no input title
saveFile = ""     # Default no saveFileName
hotFileName = ""
coldFileName = ""
keepDirectory = ""
plotFrequency = False
doScaleAve = False
doZero = False
#assume no pointing offsets
dAz = 0.
dEl = 0.
# define reference frequency for velocities (MHz)
NUH1 = 1420.405751768 # neutral hydrogen frequency (MHz)
NUOH1= 1612.231   # OH line
NUOH2= 1665.402   # OH line
NUOH3= 1667.359   # OH Line
NUOH4= 1720.530   # OH Line

# select the frequency for plotting velocities
nuRefFreq = NUH1

# define hot and cold load temperatures
thot = 285.0  # define hot and cold load temperatures
tcold = 2.7
nmedian = 0

iarg = 1
nargs = len(sys.argv)
if nargs < 3:
    print("T: Comput Tsys calibrated horn observations")
    print("Usage: T [-F <order>] [-L <velocity>] [-H <velocity>] <_seconds> <files>")
    print("Where <seconds>: Number of seconds of observations to average.")
    print("-A <hot file> <cold file> use save hot and cold load files")
    print("-B optionally plot/keep the baseline subtratcted spectra")
    print("-C optionally flag the center of the band")
    print("-D optionally print extra debugging info")
    print("-E optionally to not estimate Barycentric Velocity offset")
    print("-F <order> optionally do a polynomial baseline fit")
    print("-G <halfwidth> median filter the output vector")
    print("-H optionally set the high velocity region for baseline fit")
    print("-I optionally set Processor/Telescope Index on Plot Label")
    print("-K <directory> keep average hot and cold load calibration files")
    print("-L optionally set the low velocity region for baseline fit")
    print("-N <number> optionally set the number of spectra to plot")
    print("-M Skip writing header for .kel files")
    print("-O <AzOffset> <ElOffset> Add offsets to input Az,El")
    print("-P <dir> write PNG and PDF files instead of showing plot")
    print("-Q optionally plot intensity versus freQuency, instead of velocity")
    print("-R optionally flag known RFI lines")
    print("-S <filename> optionally set summary file name")
    print("-T <Plot Title> optionally set the plot title")
    print("-U <freqMHz> optionally set reference frequency for different line")
    print("   ie -U 1612.231, 1665.402, 1667.349, 1720.530 or 1420.40575")
    print("-V Limit ascii file putput to low and high velocity ranges")
    print("-W optionally write the calibrated Tsys files (kelvins)")
    print("-X <temp> optionally set Cold Load Temperature (Kelvins)")
    print("-Y <temp> optionally set Hot  Load Temperature (Kelvins)")
    print("-Z <file tag> optionally add tag to PDF and PNG file names")
    print("-0 optionally plot zero intensity line(s)")
    print("-MINEL optionally set the lowest elevation allowed for calibration obs (default 60d)")
    print("Observation file list must include at least one hot load file")
    print("")
    print("Glen Langston - NSF April 22, 2022")
    sys.exit()

# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if hot and cold file names are provided
    if sys.argv[iarg].upper() == '-A':
        iarg = iarg+1
        hotFileName = sys.argv[iarg]
        iarg = iarg+1
        coldFileName = sys.argv[iarg]
        print(("Calibrating with files: %s and %s" %(hotFileName, coldFileName)))
    elif sys.argv[iarg].upper() == '-B':
        doBaseline = True
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
        print("Flagging Center of Spectrum")
    elif sys.argv[iarg].upper() == '-D':
        print( 'Adding Debug Printing')
        doDebug = True
    elif sys.argv[iarg].upper() == '-E':
        print( 'Skipping Barycentric velocity correction')
        doBary = False
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
            print(("Fitting a %d-nd order polynomical baseline" % (fitOrder)))
    elif sys.argv[iarg].upper() == '-G':
        iarg = iarg+1
        nmedian = int( sys.argv[iarg])
        if nmedian < 2:
            nmedian = 2
        print("Median Filter Half Width: %d" % (nmedian))
    elif sys.argv[iarg].upper() == '-H':
        iarg = iarg+1
        maxvel = float( sys.argv[iarg])
        print(('Maximum (high) velocity for sum: %7.2f km/sec' % (maxvel)))
    elif sys.argv[iarg].upper() == '-I':
        iarg = iarg+1
        cpuIndex = float( sys.argv[iarg])
        print(('Telescope Index: %2d' % (cpuIndex)))
    elif sys.argv[iarg].upper() == '-K':
        doKeep = True
        iarg = iarg+1
        keepDirectory = sys.argv[iarg]
        keepDirectory = keepDirectory.strip()
        print('Keeping Average hot and cold files in directory: %s' % \
              (keepDirectory))
        nKeep = len(keepDirectory)
        if nKeep < 1:
            keepDirectory = "../"
        else:
            if keepDirectory[nKeep - 1] != '/':
                keepDirectory = keepDirectory + "/"
    elif sys.argv[iarg].upper() == '-L':
        iarg = iarg+1
        minvel = float( sys.argv[iarg])
        print(('Minium (low)  velocity for sum: %7.2f km/sec' % (minvel)))
    elif sys.argv[iarg].upper() == '-M':
        doWriteHeader = False
        print( 'Skipping writing file header, if writing files')
    elif sys.argv[iarg].upper() == '-N':   # if number of spectra to plot
        iarg = iarg+1
        maxPlot = int(sys.argv[iarg])
        if maxPlot < 1:
            print("Not Plotting")
        else:
            print(("Plot will have a maximum of %d spectra" % (maxPlot)))
    elif sys.argv[iarg].upper() == '-O':   # if az, el offsets
        iarg = iarg+1
        dAz = float(sys.argv[iarg])
        iarg = iarg+1
        dEl = float(sys.argv[iarg])
        print("Applying input Az,El offsets %.1f,%.1f" % (dAz, dEl))
    elif sys.argv[iarg].upper() == '-P':
        doPlotFile = True
        iarg = iarg+1
        plotFileDir = sys.argv[iarg]
    elif sys.argv[iarg].upper() == '-Q':
        plotFrequency = True
        doLimitVel = False
        print("Limiting Velocities in the plot")
    elif sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-S':   # if save file name provided
        iarg = iarg+1
        saveFile = sys.argv[iarg]
        writeTsys = True
        print( 'Saving integrated Intensities; File: %s' % ( saveFile))
    elif sys.argv[iarg].upper() == '-T':   # if plot title provided
        iarg = iarg+1
        myTitle = sys.argv[iarg]
        print('Plot Title : %s' % (myTitle))
    elif sys.argv[iarg].upper() == '-U':   # if nU ref is provided (in MHz)n
        iarg = iarg+1
        nuRefFreq = float(sys.argv[iarg])
        print(( 'Reference Frequency : %9.3f MHz' % (nuRefFreq)))
    elif sys.argv[iarg].upper() == '-V':   # limit velocity output
        doLimitVel = True
        print('Limiting Output Velocities to User Specified range')
    elif sys.argv[iarg].upper() == '-X':   # if cold load temperature
        iarg = iarg+1
        tcold = float(sys.argv[iarg])
        print('Cold Load Reference Temperature: %9.3f Kelvins' % (tcold))
    elif sys.argv[iarg].upper() == '-Y':   # if hot load
        iarg = iarg+1
        thot = float(sys.argv[iarg])
        print(( 'Hot  Load Reference Temperature: %9.3f Kelvins' % (thot)))
    elif sys.argv[iarg].upper() == '-MINEL':  # if min elevation
        iarg = iarg + 1
        lowel = float( sys.argv[iarg])
        print(( "Using elevations > %7.2f (d) for Cold load calculation" % \
                (lowel)))
    elif sys.argv[iarg].upper() == '-W':
        writeKelvin = True
    elif sys.argv[iarg].upper() == '-Z':     # label written files
        iarg = iarg+1
        fileTag = str(sys.argv[iarg])
        print(( 'File tag: %s' % (fileTag)))
    elif sys.argv[iarg].upper() == '-0':
        doZero = True
        print('Plotting zero intensity lines')
    else:
        break
    iarg = iarg + 1
    timearg = iarg
    namearg = iarg+1
# end of while not reading file names

# to create plots in cronjobs, must use a different backend
if doPlotFile:
    mpl.use('Agg')
if doBary:
    try:
        from PyAstronomy import pyasl
        BARYCENTERAVAILABLE = True
    except:
        print("!!!! Unable to import PyAstronomy !!!!")
        print("Can not compute Bary Center velocity offset")
        print("In Linux, try: ")
        print("sudo pip install PyAstronomy")
        print("or")
        print("sudo pip3 install PyAstronomy")
        BARYCENTERAVAILABLE = False
else:
    BARYCENTERAVAILABLE = False

#first argument is the averaging time in seconds
try:
    avetimesec = float(sys.argv[timearg])
except:
    print("Error: Can not parse %s as a valid average time (seconds)" % \
          sys.argv[timearg])
    sys.exit()

if nlist > nwidth and doDebug:
    print("Increasing line width array size")
    newwidth = (np.arange( 0, nlist) * 0)
    newwidth = newwidth + linewidth[nwidth-1]
    #copy over previously set values
    newwidth[0:nwidth-1] = linewidth[0:nwidth-1]
    linewidth = newwidth
    print(linewidth)
#linelist = [1400.00, 1420.0]  # RFI lines in MHz
#linewidth = [5, 5]

if plotFrequency:
    print(( "Plotting Intensity vs Frequency, Average time: %d (seconds)" % (avetimesec)))
else:
    print(( "Plotting Intensity vs Velocity, Average time: %d (seconds)" % (avetimesec)))
newObs = False
allFiles = False

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['b','r','g', 'b','r','g','b','r','g','c','m','y','c','m','y','c','m','y','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g']
nmax = len(colors)
xallmax = -9.e9
xallmin = 9.e9
ymin = 9.e9
ymax = -9.e9
yallmax = ymax
yallmin = ymin
# velocities for fitting baselines

tmin = 20.0
tmax = 999.0 # define reasoanable value limits

# prepare to remove a linear baseline

nplot = 0
nwrite = 0       # count files written out
minGlat = +90.
maxGlat = -90.
lastfreq = 0.
lastbw = 0.
lastgain = 0.
lastel = 0.
lastaz = 0.
firstdate = ""
lastdate = ""
firstaz = -1
otheraz = -1
dt = datetime.timedelta(seconds=0.)

# rest of arguments are file names
names = sys.argv[namearg:]
nFiles = len(names)
if nFiles < 2:
    print(("Not enough names to calibrate: %s" % (nFiles)))

# create the spectrum class/structure to receive spectra
rs = radioastronomy.Spectrum()
ave_cold = radioastronomy.Spectrum()
ave_hot = radioastronomy.Spectrum()

# assume only a limited range of galactic latitudes are available
# not range about +/-60.
lowGlat = 60.
# try a lower galactic latitude
lowGlat = 40.

# if a cold file provided on input
if coldFileName != "":
    ave_cold.read_spec_ast( coldFileName)
    if doScaleAve:
        ave_cold.ydataA = ave_cold.ydataA/ave_cold.count
    else:
        ave_cold.ydataA = ave_cold.ydataA/ave_cold.nave

    ncold = 1
else:
    # else must read the cold files and average
    ave_cold, minel, maxel, ncold = \
        hotcold.read_cold( names, ave_cold, lowel, lowGlat, doScaleAve)
    # if no cold obs found, first try all galactic latitudes
    if ncold < 1:
        print( "No high galactic latitude obs, trying all galactic latitudes")
        ave_cold, minel, maxel, ncold = \
            hotcold.read_cold( names, ave_cold, lowel, 0., doScaleAve)
        if ncold < 1:
            print( "No high elevation obds, using all obs > 10.")
            lowel = 10.
            ave_cold, minel, maxel, ncold = \
                hotcold.read_cold( names, ave_cold, lowel, 0., doScaleAve)
            if ncold < 1:
                print( "No high elevation data: can not calibrate")
                sys.exit()

# set number of channels for computing tSys RMS and also baseline
NRMS = 25
nData = len( ave_cold.ydataA)
# match frequency for plotting
ave_cold.refFreqHz = nuRefFreq * 1.e6
if doDebug:
    print("Ave Cold Center frequency: %7.2f MHz (%7.2f)" % \
          (ave_cold.centerFreqHz * 1.e-6, ave_cold.xdata[int(nData/2)] * 1.e-6))

# now compute ranges for baseline fitting, based on min and max vels
vel, xa0, xa, xb, xbe = \
    hotcold.velocity_indecies( ave_cold, NRMS, minvel, maxvel)

if doDebug or xa0 <=0 or xbe >= nData:
    print( 'Min Vel  %7.1f, Max Vel  %7.1f' % ( minvel, maxvel))
    print( 'Min Chan %7d, Max Chan %7d' % (xa, xb))
    print( 'Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % (minGlat, maxGlat))

if doDebug:
    print( 'Baseline fit using Vels  %7.1f to %7.1f km/s and ' % \
           ( vel[xa0], vel[xa]))
    print( 'Baseline fit using Vels  %7.1f to %7.1f km/s'  % \
           ( vel[xb], vel[xbe]))

# convert to MHz
xv = ave_cold.xdata * 1.E-6
# user input flagging uses MHZ
if flagRfi:    # if flaggin RFI
    yv = ave_cold.ydataA
    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
    ave_cold.ydataA = cv
else:
    cv = ave_cold.ydataA

if flagCenter:             # if flagging spike in center of plot
    hotcold.flagCenter( cv, nData)

if nmedian > 1:
    print("Median Filtering cold-load Obs")
    cv = tsys.medianfilter( cv, nmedian)

# initialize logging of elevation ranges
coldminel = 90
coldmaxel = 0.
minel = 200.
maxel = -200.

# read hot file or compute from the input names
if hotFileName != "":
    ave_hot.read_spec_ast(hotFileName)
    if doScaleAve:
        ave_hot.ydataA = ave_hot.ydataA/ave_hot.count
    else:
        ave_hot.ydataA = ave_hot.ydataA/ave_hot.nave
else:
    ave_hot, coldminel, coldmaxel, minGlat, maxGlat = \
        hotcold.read_hot( names, ave_hot, doScaleAve)

# convert to MHz, for flagging
xv = ave_hot.xdata * 1.E-6
ave_hot.refFreqHz = nuRefFreq * 1.e6
# do more cleanup on spectra for RFI
if flagRfi:
    yv = copy.deepcopy(ave_hot.ydataA)
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
    ave_hot.ydataA = hv
else:
    hv = copy.deepcopy(ave_hot.ydataA)

nData = len(hv)
if flagCenter:             # if flagging spike in center of plot
    hotcold.flagCenter( hv, nData)

if nmedian > 1:
    print("Median Filtering hot-load Obs")
    hv = tsys.medianfilter( hv, nmedian)

# finally compute gains
gain, gainAve, tRxMiddle, tRms, tStdA, tStdB = \
    hotcold.compute_gain( hv, cv, xa0, xa, xb, xbe, thot, tcold)

# if keeping hot and cold file names
if doKeep:
    keepExists = os.path.exists(keepDirectory)
    if keepExists:
        print("Saving average hot and cold files to existing directory")
    else:
        os.makedirs(keepDirectory)
        print("Saving average hot and cold files to new directory")
    # code for keeping hot and cold avefrage spectra
    ave_hot.ydataA = hv
    ave_cold.ydataA = cv
    if cpuIndex == 0:
        print("Update the cpu Index to distinguish telescopes (-I #)")
    hotcold.keep_hotcold( ave_hot, ave_cold, cpuIndex, keepDirectory)
    print("To Use the files for other observations, add argument:")
    print("-A %s%s %s%s " % (keepDirectory, ave_hot.filename, \
                             keepDirectory, ave_cold.filename))
print(( "Median Receiver Temp: %7.3f +/- %6.3f (%6.3f %6.3f) (K)" % \
        ( tRxMiddle, tRms, tStdA, tStdB)))

gainA = np.median(gain[xa0:xa])
gainB = np.median(gain[xb:xbe])
gainAve = 2.0/(gainA + gainB)  # Report gain in K per Count

# if plotting
if maxPlot > 0:
    fig, ax1 = plt.subplots(figsize=(10, 6))

# convert the average time in seconds into datetime format
avetime = datetime.timedelta(seconds=avetimesec)
# cound number of spectra read and averaged
nRead = 0
nave = 0

previousel = 0.

rs = radioastronomy.Spectrum()
# now read through all data and average cold sky obs
for filename in names:

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    nparts = len(parts)
    if nparts < 2:
        print("Skipping Non Astronomy file: %s" % (filename))
        continue
    aname = parts[0]
    extension = parts[1]
    nRead = nRead + 1
    if doDebug:
        print(( "%3d: %s " % (nRead, filename)))

    if nRead == nFiles:   # if this is the last file, must force output
        allFiles = True

#   create spectrum structure, empty
    rs.read_spec_ast(filename)
    if doScaleAve:
        rs.ydataA = rs.ydataA/rs.count
    else:
        rs.ydataA = rs.ydataA/rs.nave
#  apply input pointing offsets
#  assume dAz is lean east-west of the horn
    rs.utc = rs.utc + datetime.timedelta(seconds=86400.*dAz/360.)
#  stop assumming a twist in horn base
#    rs.telaz = rs.telaz + dAz
    rs.telel = rs.telel + dEl
    rs.azel2radec()    # compute ra,dec from az,el
# if not a sky observation
    if rs.telel < 0. and (not allFiles):
        continue

    # track if multiple azimuths
    if firstaz < 0:
        firstaz = rs.telaz
        otheraz = firstaz
    if firstaz != rs.telaz:
        otheraz = rs.telaz

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

    if rs.telel > maxel:
        maxel = rs.telel
    if rs.telel < minel:
        minel = rs.telel

    newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
    newFreq = (lastfreq != rs.centerFreqHz)
    newObs =  newFreq or (lastbw != rs.bandwidthHz) or \
                          (lastgain != rs.gains[0]) or newAzEl
    # if this is the last file and there was not a change in observing setup
    if allFiles and (not newObs):
        if rs.nChan != ave_spec.nChan:
            print(("File Size Error: %d != %d for file %s; Skipping ..." % \
                   (rs.nChan, ave_spec.nChan, filename)))
        else:
            ave_spec, nave, firstutc, lastutc = \
                hotcold.average_spec( ave_spec, rs, nave, firstutc, lastutc)

    # if time to average (or end of all files)
    if (dt > avetime) or newObs or allFiles:

        # not calibrating hot load observations.
        if ave_spec.telel < 0.:
            # Reset the ncold count and restart sum
            nave = 0
            continue
        ave_spec = hotcold.normalize_spec( ave_spec, firstutc, lastutc)
        xv = ave_spec.xdata * 1.E-6  # covert to MHz
        yv = ave_spec.ydataA
        if flagRfi:  # if interpolating over regions with rfi
            yv = interpolate.lines( linelist, linewidth, xv, yv)

        if nmedian > 1:
            yv = tsys.medianfilter( yv, nmedian)

        xmin = min(xv)
        xmax = max(xv)
        xallmin = min(xmin, xallmin)
        xallmax = max(xmax, xallmax)
        count = ave_spec.count
        note = ave_spec.noteA
        ncolor = min(nmax-1, nplot)

        if flagCenter:
            hotcold.flagCenter( yv, nData)

        tsky = hotcold.tsky_gain( yv, gain)

# now compute ranges for baseline fitting, based on min and max vels
        vel, xa0, xa, xb, xbe = \
            hotcold.velocity_indecies( ave_spec, NRMS, minvel, maxvel)

# get tsys from averages of ends of spectra
        tSysA = np.median(tsky[xa0:xa])
        tSysB = np.median(tsky[xb:xbe])
        tStdA = np.std(tsky[xa0:xa])
        tStdB = np.std(tsky[xb:xbe])
        cA = np.median(cv[xa0:xa])
        cB = np.median(cv[xb:xbe])
        counts = (cA+cB)/2.
        tStd = (tStdA+tStdB)/2.
        tSys = (tSysA+tSysB)/2.
        ave_spec.tSys = tSys
        ave_spec.tRx = tRxMiddle
        ave_spec.tRms = tStd
        ave_spec.tint = ave_spec.durationSec
        ave_spec.bunit = 'Kelvins'
        ave_spec.KperC = gainAve
        ave_spec.azel2radec()    # compute ra,dec from az,el and average utc

        # compute velocity correction for this direction and date
        if BARYCENTERAVAILABLE:
            corr = gf.compute_vbarycenter( ave_spec)
            if corr < minBaryVel or minBaryVel == 0.:
                minBaryVel = corr
            if corr > maxBaryVel or maxBaryVel == 0.:
                maxBaryVel = corr
            velcorr = vel + corr
        else:
            velcorr = vel
            corr = 0

        # pull out coordinates for labeling
        az = ave_spec.telaz
        el = ave_spec.telel
        gallon = ave_spec.gallon
        gallat = ave_spec.gallat
        label = 'L,L=%5.1f,%6.1f' % (gallon, gallat)
        if minGlat > gallat:
            minGlat = gallat
        if maxGlat < gallat:
            maxGlat = gallat

        nChan2 = int(ave_spec.nChan/2)
        dV = (velcorr[nChan2] - velcorr[nChan2 - 4])*.25
        if dV < 0:
            dV = - dV

        # Computes and subtracts baseline for source intensities

        # next compute the integrated intensities after baseline subtraction
        baseline = gf.fit_range( velcorr[0:nData], tsky[0:nData], \
                                 xa0, xa, xb, xbe, fitOrder)

        # remove baseline to get the source spectrum
        tSource = tsky[0:nData] - baseline[0:nData]

        # if plotting/keeping the baseline subtracted spectra, transfer to Sky
        if doBaseline:
            tsky = tSource

        # now compute integrated intensity and noise estimates
        nv = xb - xa
        if not plotFrequency:
            # create sub-arrays of intensity and velocity
            tSs = tSource[xa:xb]
            vSs = velcorr[xa:xb]
            tSourcemin = min(tSs)
            tSourcemax = max(tSs)
            # get index to maximum value; then get velocity
            iSourcemax = np.argmax(tSs)
            velSource = vSs[iSourcemax]
            # integrate over spectrum for required velocity range
            tSum = np.sum(tSs)
            tSumRms = np.std(tSs)

            # computed the integrated velocity
            tvs = tSs*vSs
            tVSum = np.sum(tvs)
            tVSumRms = np.std(tvs)

            # Integration is reported in Kelvin*Km/Sec;
            # Multiply by source velocity range
            tSumKmSec = tSum * ( maxSVel - minSVel)/float(nv)
            dTSumKmSec = tSumRms * dV * np.sqrt(float(nv))
        else:
            tSourcemin = 0.
            tSourcemax = 0.
            velSource = 0.
            tVSum = 0.
            tVSumRms = 0.
            tSumKmSec = 0.
            dTSumKmSec = 0.
            tSum = 0.

        # doppler shift in velocities channges chanels to plot
        imin, imax = gf.velocity_to_indicies( velcorr, minvel, maxvel)
        # set min and max y for plotting (only)
        nsky = len(tsky)
        if nsky < 128:
            print("Sky array too small: %d" % (nsky))
        if imin < 0:
            imin = int(nsky/8)
        if imax < imin:
            imax = int( 7*nsky/8)
        if imin > nsky:
            imin = int(nsky/8)
        if imin  > imax:
            imax = int( 7*nsky/8)

        ymin = min(tsky[imin:imax])
        ymax = max(tsky[imin:imax])

        # now compute intensity weighted velocity
        if tSum > 0.:
            tVSum = tVSum/tSum
            tVSumRms = tVSumRms*np.sqrt(float(nv))/tSum
        else:
            tVSumRms = 0.
            tVSum=0.

        # print diagnostic of integration
        if doDebug and (not plotFrequency):
            print(("Min, max Velocity   : %7.1f  to %6.1f; %d,%d" % (minvel,maxvel, imin, imax)))
            print(("Min, max Velocity I : %7.1f  to %6.1f; %d,%d" % (minSVel,maxSVel, xa, xb)))
            print(("Average Intensity(K): %7.3f +/- %6.3f K (%d)" % (tSum/nv, tSumRms/nv, nv)))
            print(("Int.  Vel.  K km/s  : %7.0f +/- %6.0f K*km/sec" % (tSumKmSec,  dTSumKmSec)))
            print(("Peak    Velocity    : %7.1f +/- %6.1f km/sec (%d)" % (velSource,  dV, iSourcemax)))
            print(("Integrated Velocity : %7.1f +/- %6.1f km/sec" % (tVSum,  tVSumRms)))

        avedatetime = ave_spec.utc.isoformat()
        datestr = avedatetime.split('T')
        atime = datestr[1]
        timeparts = atime.split('.')
        labeltime = timeparts[0]
        label = '%s, A,E: %5s,%5s, L,L: %5.1f,%6.1f' % (labeltime, az, el, gallon, gallat)
        if (minel == maxel) or (el == previousel):
            #            label = '%s L,L=%5.1f,%5.1f (%d)' % (labeltime, gallon, gallat, nave)
            el = el
        else: 
            label = '%s L,L=%5.1f,%5.1f A,E=%4.0f,%4.0f' % (labeltime, gallon, gallat, az, el)
        label = '%s L,L=%5.1f,%5.1f A,E=%4.0f,%4.0f' % (labeltime, gallon, gallat, az, el)
        if (nplot < maxPlot) or (nplot == int(nplot/20)*20):
            print(( ' Max: %9.1f  Median: %8.2f +/- %5.2f %3d %s' % \
                (tSourcemax, tSys, tStd, nave, label)))
        # if plotting frequency overwrite velocities with frequency
        if plotFrequency:
            velcorr = xv
        if int(nplot) < int(maxPlot):
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            if gallat < 7.5 and gallat > -7.5:
                plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor], \
                         label=label, lw=4)
            elif gallat < 15. and gallat > -15.:
                plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor], \
                         label=label, lw=2)
            else:
                plt.plot(velcorr, tsky, colors[ncolor], linestyle=linestyles[ncolor], \
                         label=label)
        nplot = nplot + 1

        if writeTsys:
            gf.saveTsysValues( saveFile, ave_spec, cpuIndex, tSourcemax, \
                               velSource, dV, tVSum, tVSumRms, tSumKmSec, \
                               dTSumKmSec)

        if writeKelvin:
            write_spec = copy.deepcopy( ave_spec)
            aveutc = ave_spec.utc
            outname = radioastronomy.utcToName( aveutc)
            outname = outname + ".kel"  # output in Kelvins
            write_spec.count = 1
            write_spec.nave = 1
            doComputeX = False  # x-axis already computed
            if doLimitVel:
                nOut = imax - imin
                nChan = ave_spec.nChan
                if nwrite < 1:
                    print("Limiting File Output to selected channels")
                    print("Writing Channels %d to %d (%d total)" % \
                          (imin, imax, nOut))
                nwrite = nwrite + 1
                write_spec.nChan = nOut
                # adjust x axis indicies
                i2 = int((imax+imin)/2)
                write_spec.refChan = i2 - imin
                write_spec.centerFreqHz = ave_spec.xdata[i2]
                write_spec.bandwidthHz = ave_spec.bandwidthHz*float(nOut/nChan)
                write_spec.ydataA = tsky[imin:imax]
                if plotFrequency:
                    write_spec.xdata = ave_spec.xdata[imin:imax]
                else:
                    write_spec.xdata = velcorr[imin:imax] * 1000.
                write_spec.nChan = nOut
            else:
                write_spec.ydataA = tsky
                if not plotFrequency:
                    write_spec.xdata = velcorr * 1000. # km/sec -> m/sec
            # finally write the calibrated spectrum
            if plotFrequency:
                print("Writing Intensities vs Frequency to %s%s" %
                      (keepDirectory,outname))
            else:
                print("Writing Intensities vs Velocity to %s%s" %
                      (keepDirectory,outname))
            write_spec.write_ascii_file(keepDirectory, outname, plotFrequency, \
                                       doWriteHeader, doComputeX)
        # flag restarting the sum
        nave = 0
        previousel = el

    if newObs:
        if lastfreq != rs.centerFreqHz:
            print( "Change: LastFreq: %8.1f  New  %8.1f (MHz)" % \
                    (lastfreq/1e6, rs.centerFreqHz/1e6))
            lastfreq = rs.centerFreqHz
        if lastbw != rs.bandwidthHz:
            print( "Change: LastBandwidth: %5.1f  New: %5.1f (MHz)" % \
                    (lastbw/1e6, rs.bandwidthHz/1e6))
            lastbw = rs.bandwidthHz
        if lastgain != rs.gains[0]:
            print( "Change: LastGain: %4.1f  New: %4.1f (db) " %\
                    (lastgain, rs.gains[0]))
            lastgain = rs.gains[0]
        if newAzEl:
            print( "Change: LastAzEl: %6.1f,%6.1f  New: %6.1f,%6.1f (deg) " %\
                    (lastaz,lastel, rs.telaz, rs.telel))
            lastaz = rs.telaz
            lastel = rs.telel

    # end if a cold file
    if allFiles:
        break

    if rs.telel > 0:
    # Average in most recently read spectrum
        if rs.nChan != ave_spec.nChan:
            print(("File Size Error: %d != %d for file %s; Skipping ..." % (rs.nChan, ave_spec.nChan, filename)))
        else:
            ave_spec, nave, firstutc, lastutc = \
                hotcold.average_spec( ave_spec, rs, nave, firstutc, lastutc)

    #end for all files to sum
# end of all files to read

print(("Min, Max Galactic Latitude: %7.1f,%7.1f" % (minGlat, maxGlat)))
print(("Min, Max Elevation:         %7.1f,%7.1f" % (minel, maxel)))

if doDebug:
    print(( 'Number of remaining observations not plotted: %d' % ( ncold)))

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
    myTitle = myTitle + ("I:%d " % (cpuIndex))

if firstaz == otheraz:
    myTitle = myTitle + "Az = %6.1f, " % (firstaz)
else:
    myTitle = myTitle + "Az = %6.1f to %6.1f, " % (firstaz, otheraz)

if minel == maxel:
    myTitle = myTitle + " El=%6.1f" % (minel)
else:
    myTitle = myTitle + " El=%6.1f to %6.1f" % (minel, maxel)

if (maxPlot < 1) or (nplot < 1):
    print("No Plots, exiting")
    sys.exit()
#fig.canvas.set_window_title(myTitle)
#for tick in ax1.xaxis.get_major_ticks():
#    tick.label.set_fontsize(14)
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small')
    #tick.label.set_rotation('vertical')
#for tick in ax1.yaxis.get_major_ticks():
#    tick.label.set_fontsize(14)
#plt.xlim(xallmin,xallmax)  # normal plot range
#plt.xlim(-150., 300.)

plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

# select the relevant range of velocities (km/sec) for plotting
#plt.xlim(-250., 250.)
if not plotFrequency:
    plt.xlim(minvel, maxvel)

if doZero:
    ax1.axhline( linewidth=1, linestyle=':', color='k')

dy = yallmax - yallmin
if dy < 8:
    dy = 8

if doDebug:
    print(( "Y min, max: %8.3f, %8.3f" % (yallmin, yallmax)))

# set the y scale to go above and below all data
plt.ylim((yallmin-(dy/8.)), (yallmax+(dy/4.)))
plt.title(myTitle, fontsize=16)
if plotFrequency:
    plt.xlabel('Frequency (MHz)', fontsize=16)
else:
    if minBaryVel != 0. or maxBaryVel != 0.:
        if minBaryVel != maxBaryVel:
            plt.xlabel(('Velocity (km/sec; V_corr: %.1f to %.1f)' % \
                        (minBaryVel, maxBaryVel)), fontsize=16)
        else:
            plt.xlabel(('Velocity (km/sec; V_corr:%.1f)' % corr), fontsize=16)
    else:
        plt.xlabel('Velocity (km/sec)', fontsize=16)

plt.ylabel('Intensity (Kelvins)', fontsize=16)
plt.legend(loc='upper right')
# if writing files
if doPlotFile:
    if fileTag == "":
        fileTag = "T-" + firstdate
    outpng = plotFileDir + fileTag + ".png"
    plt.savefig(outpng,bbox_inches='tight')
    outpdf = plotFileDir + fileTag + ".pdf"
    plt.savefig(outpdf,bbox_inches='tight')
    print(( "Wrote files %s and %s" % (outpng, outpdf)))
else:
    # else show the plots
    plt.show()
