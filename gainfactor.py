"""
functions to compute gain factors for different processor indicies, dates
and elevations.  This module also finds the local times of galactic plane
Crossings.
"""

# HISTORY
# 24Oct15 GIL improve check for blank lines and end of file
# 24Mar02 GIL enable returning telescope longitude and latitude
# 23Oct26 GIL add another digit to tSys, tRms fix label
# 23Oct10 GIL add RA,Dec and Galactic Longitude and Latitude to log file
# 23Aug16 GIL add velocity calculation
# 22May04 GIL fix reading LONLAT etc
# 22May02 GIL revise save file; add telescope location
# 21AUG20 GIL enable passing of the debug flag
# 21JAN07 GIL add help, update for changes
# 20APR30 GIL add pointing offset model
# 20APR29 GIL add read savefile function
# 20APR22 GIL clean up save labels
# 20APR21 GIL add intensity weighted velocity error
# 20APR20 GIL add intensity weighted velocity moment
# 20APR17 GIL add peak source in velocity range to log
# 20APR16 GIL add reading and writing of gain factor logs
# 20APR06 GIL initial version based on holdcold.py
#
import os.path
import datetime
import numpy as np

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

# define a few parameters
NMEAS = 10
dates = np.zeros(NMEAS)
els = np.zeros(NMEAS)
azs = np.zeros(NMEAS)
nMeas = np.zeros(NMEAS)
CLIGHT = 299792458.0  # speed of light in m/sec

nObs = 1
nProcessors = 8
nMeasurements = 3
measurements = np.zeros((nObs, nProcessors, nMeasurements))
# now fill in measuremnts for this observations
iObs = 0
azs[iObs] = 180.0
els[iObs] = 90.0
nMeas[iObs] = 3
measurements[
    iObs,
    2,
] = (26.3, 26.0, 25.3)
measurements[
    iObs,
    3,
] = (49.6, 47.7, 44.9)
measurements[
    iObs,
    4,
] = (47.2, 48.2, 48.4)
measurements[
    iObs,
    5,
] = (40.0, 39.6, 38.1)
gainFactors = np.zeros((nObs, nProcessors))
normalized = False
firstRun = True

# save global position variables
lastAz = -100.0
lastEl = -100.0
lastId = 0

# telescope azimuth and elevation (degrees)
telEl = 0.0
telAz = 0.0
# telescope location (degrees)
telLon = -79.831111111
telLat = 38.422500000
telAlt = 820.000

def fit_baseline(xs, ys, imin, imax, nchan, fitOrder, doDebug=False):
    """
    fit baseline does a polynomical fit over channels in a select range
    The baseline is returned.
    Inputs:
    xs     x axis values
    ys     y axis values
    imin   index to center location of iminimum side fit
    imax   index to center location of imaximum side fit
    nchan
    """

    xfit = np.concatenate(
        (xs[imin - nchan : imin + nchan], xs[imax - nchan : imax + nchan])
    )
    yfit = np.concatenate(
        (ys[imin - nchan : imin + nchan], ys[imax - nchan : imax + nchan])
    )
    # calculate polynomial (0=constant, 1=linear, 2=2nd order, 3=3rd order
    z = np.polyfit(xfit, yfit, fitOrder)
    if doDebug:
        print("2nd order Fit Coefficients: %s" % (z))
    f = np.poly1d(z)  # f is a function that is used to compute y values

    # calculate y's from the xs
    yout = f(xs)

    return yout


def fit_range(xs, ys, xa0, xa, xb, xbe, fitOrder, doDebug=False):
    """
    fit baseline does a polynomical fit over channels in a select range
    The baseline is returned.
    Inputs:
    xs     x axis values
    ys     y axis values
    xa0, xa  one
    xb, xbe  2nd range of x,y channels to fit (maximum side)
    nchan
    """

    xfit = np.concatenate((xs[xa0:xa], xs[xb:xbe]))
    yfit = np.concatenate((ys[xa0:xa], ys[xb:xbe]))
    # calculate polynomial (0=constant, 1=linear, 2=2nd order, 3=3rd order
    z = np.polyfit(xfit, yfit, fitOrder)
    if doDebug:
        print("2nd order Fit Coefficients: %s" % (z))
    f = np.poly1d(z)  # f is a function that is used to compute y values

    # calculate y's from the xs
    yout = f(xs)

    return yout


def compute_vbarycenter(spectrum, doDebug=False):
    """
    Compute the velocity correction to Barycentric for this date and direction
    """
    global firstRun

    if baryCenterAvailable:
        telLon = spectrum.tellon
        telLat = spectrum.tellat
        altitude = spectrum.telelev
        ra2000 = spectrum.ra
        dec2000 = spectrum.dec

        # need to convert date time to Julian date
        jd = jdutil.datetime_to_jd(spectrum.utc)

        # Calculate barycentric correction (debug=True show
        # various intermediate results)
        corr, hjd = pyasl.helcorr(
            telLon, telLat, altitude, ra2000, dec2000, jd, debug=False
        )
        if doDebug:
            print("Barycentric correction [km/s]: %8.3f" % (corr))
        elif firstRun:
            print("1st Barycentric correction [km/s]: %8.3f" % (corr))
        firstRun = False
    else:
        corr = 0.0
    return corr


def velocity_to_indicies(vel, minvel, maxvel):
    """
    Function to compute indices from velocity array and target velocities
    """

    nData = len(vel)
    iref = int(nData / 2)

    vref = vel[iref]
    dv = (vel[iref + 2] - vel[iref - 2]) / 4.0
    imin = int(((minvel - vref) / dv) + iref)
    imax = int(((maxvel - vref) / dv) + iref) + 1

    if imax < imin:  # swap indices if frequency opposit velocity
        temp = imin
        imin = imax
        imax = temp

    if imin < 0:
        print("Imin Error computing baseline: ", imin)
        imin = 0
    if imin >= nData:
        print("Imin Error computing baseline: ", imin)
        imin = nData - 1

    if imax < 0:
        print("Imax Error computing baseline: ", imax)
        imax = 0
    if imax >= nData:
        print("Imax Error computing baseline: ", imax)
        imax = nData - 1

    # some channels will be selected using x[imin:imax] so,
    # now make sure range increases
    if imin > imax:
        temp = imin
        imin = imax
        imax = temp
    return imin, imax


def velocity(freqHz, refFreqHz):
    """
    Function to compute velocities (m/s) from input frequencies
    The input frequencies are defined as Hz, but since only the ratios are used,
    Any frequency units are acceptable.
    """

    nData = len(freqHz)
    refFreqHz = float(refFreqHz)
    velmps = np.zeros(nData) + refFreqHz  # velocity in meters/second

    # doppler shift is defined as higher frequency is negative velocity
    # negative velocity means moving towards us
    for iii in range(nData):
        velmps[iii] = (CLIGHT / refFreqHz) * (velmps[iii] - freqHz[iii])

    return velmps  # end of velocity


def saveTsysValues(
    saveFile,
    cSpec,
    cpuIndex,
    tSourcemax,
    velSource,
    dV,
    tVSum,
    tVsumRms,
    tSumKmSec,
    dTSumKmSec,
):
    """
    saveTsysValues - saves the calibration values for this calibrated observation
    where
    saveFile - Name of file used for recording observations
    spectrum - Tsys Calibrated Spectrum to be saved.
    cpuIndex - index number of horn/cpu used to make the measurements
    tSourcemax - peak of intensity in identified velocity range

    The Tsys values are logged in 3 parts,
    1. Tsys value for the observation.
    2. The Counts per Kevin Value for the observation
    3. The normalization factor to put all horns on the same scale.
    The first 2 values are determined by hot and cold load observations.
    The hot ground and empty sky.
    The normaization factor (3) is measured by comparison of simultaneous astronomical observations
    of the same region of the sky.   This can be done by a 24 hour drift scan with all horns.
    """
    global lastAz, lastEl
    global telLon, telLat, telAlt

    autc = str(cSpec.utc)  # ascii version of UTC
    parts = autc.split(" ")
    nparts = len(parts)
    date = parts[0]
    nd = len(date)
    date = date[2:nd]  # trim off 2019 to 19
    time = parts[1]
    time = time.replace("_", ":")  # put time back in normal hh:mm:ss format
    parts = time.split(".")  # trim off seconds part of time
    time = parts[0]
    if saveFile == "":
        saveFile = "../" + date + "-" + str(cpuIndex), strip() + ".sav"

    # determine if the file exists already
    oldFile = os.path.isfile(saveFile)

    f = open(saveFile, "a+")

    # update coordinates and position every time an az,el change
    if lastAz != cSpec.telaz or lastEl != cSpec.telel:
        f.write(
            "#LONLAT %02d %15.9f %15.9f %9.3f \r\n"
            % (cpuIndex, cSpec.tellon, cSpec.tellat, cSpec.telelev)
        )
        f.write("#AZEL   %02d %7.2f %7.2f \r\n" % (cpuIndex, cSpec.telaz, cSpec.telel))
        lastAz = cSpec.telaz
        lastEl = cSpec.telel
        lastid = cpuIndex

    # if a new file, then need to add the header
    if not oldFile:
        f.write(
            "#  Date    Time   Tel  Az     El     Tsys    Trx    Trms     Time  K/Count    Peak   Peak  Vel.     Sum Vel.   "
        )
        f.write("Sum Intensity   Scale     RA      Dec   G Lon   GLat \r\n")
        f.write(
            "#                  #   (d)    (d)     (K)     (K)    (K)     (s)              (K)     (km/s) +/-    (km/s) +/-  "
        )
        f.write(" (K km/s) +/-   Factor")
        f.write("(d)     (d)    (d)   (d) \r\n")

    #          1 2   3   4     5     6    7      8       9    10     11    11    12     13    14      15    16     17
    #         Date   cpu  az    el  tSys  tRx   tRms   tint   K/C   tPeak, vel    dv   VSum   Vrms,  KInt  dKInt factor
    f.write(
        "%s %s %2d %6.1f %6.1f %7.2f %7.3f %6.3f %7.0f %7.1f %7.3f %7.1f %5.1f %7.1f %5.1f %7.0f %7.0f %7.3f %7.2f %7.2f %7.2f %7.2f\r\n"
        % (
            date,
            time,
            cpuIndex,
            cSpec.telaz,
            cSpec.telel,
            cSpec.tSys,
            cSpec.tRx,
            cSpec.tRms,
            cSpec.tint,
            cSpec.KperC,
            tSourcemax,
            velSource,
            dV,
            tVSum,
            tVsumRms,
            tSumKmSec,
            dTSumKmSec,
            cSpec.gainFactor,
            cSpec.ra,
            cSpec.dec,
            cSpec.gallon,
            cSpec.gallat,
        )
    )
    f.close()
    # end of saveTsysValues()


def readSummaryTelescopeLocation(f):
    """
        readSummaryTelescopeLocation steps through the initial lines of a
        'T' program Summary file and extracts the telescope geographic
        longitude and latitude.   Examples:
    #LONLAT 07   -79.831111111    38.422500000   820.000
    #AZEL   07  180.00   83.04
    """
    # link to global values
    global telLon, telLat, telAlt, telAz, telEl
    count = 0
    
    # read through the file header info
    while True:
        try:
            aline = f.readline()
        except:
            if count < 5:
                print("Unable to read Summary File")
            break
        aline = aline.strip()
        alen = len(aline)
        count = count + 1
        
        # accept blank lines or lines starting with #
        if alen > 0:
            if aline[0] != "#":
                break
            # search for keywords
            if "LONLAT" in aline[1:8]:
                values = aline[7:]
                parts = values.split()
                if len(parts) == 4:
                    cpuIndex = int(parts[0])
                    telLon = float(parts[1])
                    telLat = float(parts[2])
                    telAlt = float(parts[3])
                else:
                    print("Invalid Telescope Lon lat values:")
                    print(aline)
                    print(parts)
            if "AZEL" in aline[1:8]:
                values = aline[6:]
                parts = values.split()
                if len(parts) == 3:
                    cpuIndex = int(parts[0])
                    telAz = float(parts[1])
                    telEl = float(parts[2])
                else:
                    print("Invalid Telescope Az,El values:")
                    print(aline)
                    print(parts)

    # end of readSummaryTelescopeLocation()
    # return the last data line as well as the parsed values
    return telLon, telLat, telAlt, telAz, telEl, cpuIndex, aline


def readSaveValues(aline):
    """
    readSaveValues - saves the calibration values for this calibrated observation
    where
    f - Open file
    outputs:

    The Tsys values are logged in 3 parts,
    1. Tsys value for the observation.
    2. The Counts per Kevin Value for the observation
    3. The normalization factor to put all horns on the same scale.
    The first 2 values are determined by hot and cold load observations.
    That the hot ground and empty sky.
    The normaization factor (3) is measured by comparison of observations
    of the same region of the sky.   
    This can be done by a 24 hour drift scan with all horns.
    """

    #    print("Entering Read Save Values")
    global lastAz, lastEl

    date = ""
    time = ""
    id = int(0)
    telaz = float(0.0)
    telel = float(90.0)
    tSys = float(100.0)
    tRx = float(50.0)
    tRms = float(1.0)
    tint = float(1.0)
    KperC = float(100.0)
    tSourcemax = float(10.0)
    velSource = float(0.0)
    dV = float(1.0)
    tVSum = float(0.0)
    tVsumRms = float(1.0)
    tSumKmSec = float(1000.0)
    dTSumKmSec = float(10.0)
    gainFactor = float(1.0)
    ra = float(0.0)
    dec = float(0.0)
    gallon = float(0.0)
    gallat = float(0.0)

    # if a new file, then need to add the header
    #   f.write( "#  Date    Time   Tel  Az     El     Tsys    Trx    Trms     Time  K/Count    Peak    Peak  Vel.     Sum Vel.   ")
    #   f.write( "Sum Intensity   Scale \r\n")
    #   f.write( "#                  #   (d)    (d)     (K)     (K)    (K)     (s)              (K)     (km/s) +/-    (km/s) +/-  ")
    #   f.write( " (K km/s) +/-   Factor\r\n")

    #          1 2   3   4     5     6    7      8       9    10     11    11    12     13    14      15    16     17
    #         Date   cpu  az    el  tSys  tRx   tRms   tint   K/C   tPeak, vel    dv   VSum   Vrms,  KInt  dKInt factor
    # if hear, not a comment
    parts = aline.split()
    nparts = len(parts)
    # if at end of file no values, there are no parts
    date = ""
    if nparts < 17:
        return (
            date,
            time,
            id,
            telaz,
            telel,
            tSys,
            tRx,
            tRms,
            tint,
            KperC,
            tSourcemax,
            velSource,
            dV,
            tVSum,
            tVsumRms,
            tSumKmSec,
            dTSumKmSec,
            gainFactor,
            ra,
            dec,
            gallon,
            gallat,
        )

    date = parts[0]
    time = parts[1]
    cpuIndex = int(parts[2])
    telaz = float(parts[3])
    telel = float(parts[4])
    tSys = float(parts[5])
    tRx = float(parts[6])
    tRms = float(parts[7])
    tint = float(parts[8])
    Kperc = float(parts[9])
    tSourcemax = float(parts[10])
    velSource = float(parts[11])
    dV = float(parts[12])
    tVSum = float(parts[13])
    tVSumRms = float(parts[14])
    tSumKmSec = float(parts[15])
    dTSumKmSec = float(parts[16])
    gainFactor = float(parts[17])
    ra = float(parts[18])
    dec = float(parts[19])
    gallon = float(parts[20])
    gallat = float(parts[21])

    #        sd( "%s %s %2d %6.1f %6.1f %7.2f %7.2f %6.2f %7.0f %6.1f %7.3f %7.1f %5.1f %7.1f %5.1f %7.0f %7.0f %7.3f %7.2f %7.2f %7.2f %7.2f \r\n" %
    #             (date, time, cpuIndex, telaz, telel, tSys, tRx, tRms, tint, KperC,
    #              tSourcemax, velSource, dV, tVSum, tVsumRms, tSumKmSec, dTSumKmSec, gainFactor))

    return (
        date,
        time,
        cpuIndex,
        telaz,
        telel,
        tSys,
        tRx,
        tRms,
        tint,
        KperC,
        tSourcemax,
        velSource,
        dV,
        tVSum,
        tVsumRms,
        tSumKmSec,
        dTSumKmSec,
        gainFactor,
        ra,
        dec,
        gallon,
        gallat,
    )
    # end of readSaveValues()


def readAllValues(filename):
    """
    readAllValues - reads one file of saved intensity values and return arrays
    input is a file name for a horn radio telescope summary file
    (Usually created by the 'T' program)
    Outputs are:
    date     - ascii string date of the summary data (ie 20-05-03)
    utcs     - UTCs of measurements
    secs     - array of seconds past midnight of the first date
    tSums    - array of integrated intensity for the provided spectrum
    dTs      - array of uncertainty estimates for the integrated intensities
    gallons  - array of galactic longitudes
    gallats  - array of galactic latitudes
    """
    f = open(filename, "r")

    #    timefmt = "%Y-%m-%d %H:%M:%S"
    timefmt = "%Y-%m-%d %H:%M:%S"
    global lastId, lastAz, lastEl
    global telLon, telLat, telAlt, telAz, telEl, cpuIndex

    count = 0
    utcs = []
    secs = []
    tSums = []
    dTs = []
    tMaxs = []
    azs = []
    els = []
    ras = []
    decs = []
    gallons = []
    gallats = []
    firstdate = ""

    telLon, telLat, telAlt, telAz, telEl, cpuIndex, aline = \
        readSummaryTelescopeLocation(f)

    print("Read Telescope Info, Id: %d" % (cpuIndex))

    # read though all values in the summary file
    nblank = 0
    lastEl = -99.
    while True:
        if len(aline) < 0:
            aline = f.readline()
            aline.strip()
        # allow a few blank lines in file.  
        if len(aline) == 0:
            nblank + 1
            # if too many blank lines, must be end of file.
            if nblank > 10:
                break # exit loop
            continue  # else keep looking for next data line
            
        date, time, id, telaz, telel, tSys, tRx, tRms, tint, KperC, \
            tSourcemax, velSource, dV, tVSum, tVsumRms, tSumKmSec, \
            dTSumKmSec, gainFactor, ra, dec, gallon, gallat = \
                    readSaveValues(aline)
        if telel != lastEl:
            print("readvalues: New El %7.2f -> %7.2f" % (lastEl, telel))
            lastEl = telel
        nblank = 0
        
        # if date is blank, then end of file
        if date == "":
            break
        if firstdate == "":
            firstdate = date
            lastId = cpuIndex
            lastAz = telaz
            lastEl = telel
        if count == 0:
            utcmidnight = datetime.datetime.strptime("20" + date + " 00:00:00", timefmt)

        utc = datetime.datetime.strptime("20" + date + " " + time, timefmt)
        dUtc = utc - utcmidnight
        utcs.append(utc)
        secs.append(dUtc.total_seconds())
        tSums.append(tSumKmSec)
        dTs.append(dTSumKmSec)
        tMaxs.append(tSourcemax)
        gallons.append(gallon)
        gallats.append(gallat)
        azs.append(telaz)
        els.append(telel)
        ras.append(ra)
        decs.append(dec)
        if count > 0:
            if telel != lastEl:
                print("Telescope Elevation Changed %8.1f => %8.1f" % \
                      (lastEl, telel))
        lastEl = telel

        if count < 3:
            print("%d: %s %9.3f %7.3f" % (count, utc, tSumKmSec, dTSumKmSec))
        count = count + 1
        aline = f.readline()
        aline.strip()

    secs = np.asarray(secs)
    tSums = np.asarray(tSums)
    dTs = np.asarray(dTs)
    tMaxs = np.asarray(tMaxs)
    azs = np.asarray(azs)
    els = np.asarray(els)
    ras = np.asarray(ras)
    decs = np.asarray(decs)
    gallons = np.asarray(gallons)
    gallats = np.asarray(gallats)
    # ascii string, datetime, seconds, Temperature sum, temp uncertainty
    #    return firstdate, utcs, secs, tSums, dTs, azs, els, ras, decs, gallons, gallats
    return firstdate, utcs, secs, tMaxs, dTs, azs, els, ras, decs, gallons, gallats, cpuIndex


def readTsysValues(saveFile, utc, cpuIndex, az, el):
    """
    readTsysValues - read the calibration values for this calibrated observation
    where
    saveFile - Name of file used for recording observations
    spectrum - Tsys Calibrated Spectrum to be saved.
    cpuIndex - index number of horn/cpu used to make the measurements
    """

    autc = str(calibratedSpectrum.utc)  # ascii version of UTC
    parts = autc.split(" ")
    date = parts[0]
    nd = len(date)
    date = date[2:nd]  # trim off 2019 to 19
    time = parts[1]
    time = time.replace("_", ":")  # put time back in normal hh:mm:ss format
    parts = time.split(".")  # trim off seconds part of time
    time = parts[0]
    az = rs.telaz
    el = rs.telel
    if saveFile == "":
        saveFile = date + ".sav"

    f = open(saveFile, "w+")
    f.write("%s %s %2d %7.2f %7.2f\r\n" % (date, time, cpuIndex, az, el))
    f.close()


def lonlatelev():
    """
    latlonel() returns the last read telescope geographic location
    """
    global telLon, telLat, telAlt
    return telLon, telLat, telAlt


def lastAzel():
    """
    latlonel() returns the last read telescope geographic location
    """
    global lastId, lastAz, lastEl
    return lastId, lastAz, lastEl


def normalizemeasures(iObs):
    """
    normalizemeasures normalizes one set of gain measurements for all processors
    iObs: Is the observation set to normalize
    """
    global nMeas
    gainSum = 0.0
    if iObs < 0 or iObs >= nObs:
        print("Invalid measurement index: %d)" % (iObs))
        return 0
    n = int(nMeas[iObs])
    count = 0  # count total number of measurements
    for jP in range(nProcessors):
        for iN in range(n):
            if measurements[iObs, jP, iN] > 0:
                gainSum = gainSum + measurements[iObs, jP, iN]
                count = count + 1
    # if any measurements, compute average for all measurements
    if count > 0:
        gainAve = gainSum / float(count)
    else:
        print("No Measurements for observation Index %d" % (iObs))

    # now for each processor, compute the factor that normaizes the average
    for jP in range(nProcessors):
        gainSum = 0.0  # now computing only average fro this processor
        count = 0
        for iN in range(n):  # averages et of measuremnts
            if measurements[iObs, jP, iN] > 0:
                gainSum = gainSum + measurements[iObs, jP, iN]
                count = count + 1
        if count > 0:
            gainSum = gainSum / float(count)
            gainFactors[iObs, jP] = gainAve / gainSum
            print(
                "Obs %d, Processor %d factor: %7.2f" % (iObs, jP, gainFactors[iObs, jP])
            )
        else:
            gainFactors[iObs, jP] = 1.0  # default is unity gain

    return n


def compute_gain_factor(pIndex, aveutc, az, el):
    """
    Compute a processor based gain factor for inputs
    Inputs:
    pIndex          Index to processor based measurements
    aveutc          Date of observation for which the factors are measured
    az              Azimuth (degrees) for which the factor is computed
    el              Elevation (degrees) for which the factor is computed
    """
    global normalized

    # temporarily just use one set of values
    if not normalized:
        for iObs in range(nObs):
            normalizemeasures(iObs)
        normalized = True

    iObs = 0
    gainFactor = gainFactors[iObs, pIndex]
    return gainFactor


# end of def compute_gain_factor()


def gainfactor(aveutc, pIndex, az, el):
    """
    gain_factor returnes kelvins/count for a date-utc, processor index, azimuth and elevation)
    Inputs:
    aveutc          Date of observation utc in modified julian date
    pIndex          Index to processor based measurements
    azimuth         Azimuth (degrees)
    el              Elevation (degrees)
    Output:
    gain_factor     in Kelvins/Count
    """
    global normalized

    # temporarily just use one set of values
    if not normalized:
        for iObs in range(nObs):
            normalizemeasures(iObs)
        normalized = True

    # gain factor does not yet include azimuth or elevation corrections
    gainFactor = gainFactors[iObs, pIndex]
    return gainFactor


# end of def gainfactor()

def listSave(savefile):
    """
    listSave lists all entries in a provided savefile
    Inputs
    savefile - name of the file containing the saved observating summaries
    """

    f = open(savefile, "r")
    date = "Uknown"
    count = 0
    while date != "":
        aline = f.readline()
        date, time, id, telaz, telel, tSys, tRx, tRms, tint, KperC, \
        tSourcemax, velSource, dV, tVSum, tVsumRms, tSumKmSec, dTSumKmSec, \
        gainFactor, ra, dec, gallon, gallat = readSaveValues(aline)
        dlen = len(date)
        if dlen < 1:
            break
        if (date[0] != "#") and (count < 3):
            print(
                "%3d: %s %d %8.2f %7.2f %5.2f %7.0f"
                % (count, date, id, tSys, tRx, tRms, tint)
            )
        count = count + 1
        
    f.close()
    # end of read
    return count
