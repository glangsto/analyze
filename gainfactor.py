#Python function to compute gain factor for different processor indicies, dates and elevations
#HISTORY
#20APR17 GIL add peak source in velocity range to log
#20APR16 GIL add reading and writing of gain factor logs
#20APR06 GIL initial version based on holdcold.py
#
import sys
import datetime
import numpy as np
import radioastronomy

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
nObs = 1
nProcessors = 8
nMeasurements = 3
measurements = np.zeros((nObs, nProcessors, nMeasurements))
# now fill in measuremnts for this observations
iObs = 0
azs[iObs] = 180.
els[iObs] = 90.
nMeas[iObs] = 3
measurements[iObs, 2,] = ( 26.3, 26.0, 25.3)
measurements[iObs, 3,] = ( 49.6, 47.7, 44.9)
measurements[iObs, 4,] = ( 47.2, 48.2, 48.4)
measurements[iObs, 5,] = ( 40.0, 39.6, 38.1)
gainFactors = np.zeros((nObs, nProcessors))
normalized = False
doDebug = False
firstRun = True

def fit_baseline( xs, ys, imin, imax, nchan, fitOrder):
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
    
    xfit = np.concatenate( (xs[imin-nchan:imin+nchan],xs[imax-nchan:imax+nchan]))
    yfit = np.concatenate( (ys[imin-nchan:imin+nchan],ys[imax-nchan:imax+nchan]))
# calculate polynomial (0=constant, 1=linear, 2=2nd order, 3=3rd order
    z = np.polyfit(xfit, yfit, fitOrder)
    if doDebug:
        print("2nd order Fit Coefficients: %s" % (z))
    f = np.poly1d(z)   # f is a function that is used to compute y values

# calculate y's from the xs
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

def saveTsysValues( saveFile, cSpec, cpuIndex, tSourcemax, velSource, integratedK, rmsIntegratedK):
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
    That the hot ground and empty sky.
    The normaization factor (3) is measured by comparison of simultaneous astronomical observations 
    of the same region of the sky.   This can be done by a 24 hour drift scan with all horns.
    """

    autc = str(cSpec.utc) # ascii version of UTC
    parts = autc.split(' ')
    date = parts[0]
    nd = len(date)
    date = date[2:nd]                  # trim off 2019 to 19
    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format
    parts = time.split('.')  # trim off seconds part of time
    time = parts[0]
    if saveFile == "":
        saveFile = "../" + date + ".sav"

    f = open(saveFile, "a+")
    #          1 2   3   4     5     6    7      8      9     10       11    12    13
    #         Date   cpu  az    el  tSys  tRx   tRms tSource vSource tint   K/C   factor
    f.write( "%s %s %2d %6.1f %6.1f %7.2f %7.2f %6.2f %7.3f %6.1f %6.0f, %7.2f %9.1f %7.1f %6.3f\r\n" % 
             (date, time, cpuIndex, cSpec.telaz, cSpec.telel, cSpec.tSys, cSpec.tRx, cSpec.tRms, tSourcemax, velSource,
              cSpec.tint, cSpec.KperC, integratedK, rmsIntegratedK, cSpec.gainFactor))
    f.close()
    # end of saveTsysValues()

def readTsysValues( saveFile, utc, cpuIndex, az, el):
    """ 
    readTsysValues - read the calibration values for this calibrated observation
    where
    saveFile - Name of file used for recording observations
    spectrum - Tsys Calibrated Spectrum to be saved.
    cpuIndex - index number of horn/cpu used to make the measurements
    """

    autc = str(calibratedSpectrum.utc) # ascii version of UTC
    parts = autc.split(' ')
    date = parts[0]
    nd = len(date)
    date = date[2:nd]                  # trim off 2019 to 19
    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format
    parts = time.split('.')  # trim off seconds part of time
    time = parts[0]
    az = rs.telaz 
    el = rs.telel 
    if saveFile == "":
        saveFile = date + ".sav"

    f = open(saveFile, "w+")
    f.write( "%s %s %2d %7.2f %7.2f\r\n" % (date, time, cpuIndex, az, el))
    f.close()

def normalizemeasures( iObs):
    """
    normalizemeasures normalizes one set of gain measurements for all processors
    iObs: Is the observation set to normalize
    """
    global nMeas
    gainSum = 0.
    if (iObs < 0 or iObs >= nObs):
        print ("Invalid measurement index: %d)" % (iObs))
        return n
    n = int(nMeas[iObs])
    count = 0       # count total number of measurements
    for jP in range (nProcessors):
        for iN in range( n):
            if measurements[iObs, jP, iN] > 0:
                gainSum = gainSum + measurements[iObs, jP, iN]
                count = count + 1
    # if any measurements, compute average for all measurements
    if count > 0:
        gainAve = gainSum/float(count)  
    else:
        print ("No Measurements for observation Index %d" % (iObs))

    # now for each processor, compute the factor that normaizes the average
    for jP in range (nProcessors):
        gainSum = 0.                 # now computing only average fro this processor
        count = 0
        for iN in range( n):         # averages et of measuremnts
            if measurements[iObs, jP, iN] > 0:
                gainSum = gainSum + measurements[iObs, jP, iN]
                count = count + 1
        if count > 0:
            gainSum = gainSum/float(count)  
            gainFactors[ iObs, jP] = gainAve/gainSum
            print ("Obs %d, Processor %d factor: %7.2f" % (iObs, jP, gainFactors[iObs, jP]))
        else:        
            gainFactors[ iObs, jP] = 1.0      # default is unity gain

    return

def compute_gain_factor( pIndex, aveutc, az, el):
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
        for iObs in range (nObs):
            normalizemeasures( iObs)
        normalized = True

    iObs = 0
    gainFactor = gainFactors[iObs, pIndex]
    return gainFactor
# end of def compute_gain_factor()

def gainfactor( aveutc, pIndex, az, el):
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
        for iObs in range (nObs):
            normalizemeasures( iObs)
        normalized = True

    # gain factor does not yet include azimuth or elevation corrections
    gainFactor = gainFactors[iObs, pIndex]
    return gainFactor

# end of def gainfactor()

