"""
Module hotcold.py contains procedures for averaging files read
in from lists of names.
Glen Langston National Science Foun dation
"""
#Python Script to plot calibrated  NSF spectral integration data.
#plot the raw data from the observation
#HISTORY
#22APR22 GIL initial version based on t.py
#
import copy
import sys
import numpy as np
import radioastronomy
import gainfactor as gf

# define a small number
EPSILON = 1.E-10

C = 299792.458  # (Speed of light  km/sec)

# prepare to compute weighted averages
def average_spec( ave_spec, in_spec, nave, firstutc, lastutc):
    """
    Compute weighted average spectra
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

    # if restarting the sum
    if nave == 0:
        ave_spec = copy.deepcopy(in_spec)  # initial spectrum is one just read
        firstutc = in_spec.utc
        lastutc = in_spec.utc
        nave = 1
        # sums are weighted by durations
        ave_spec.ydataA = (in_spec.ydataA * in_spec.durationSec)  # duration scaling
        # keep track of observing time for weighted sum
        ave_spec.durationSec = in_spec.durationSec
    else: # else not enough time yet, average ave_spec data
        if in_spec.utc < firstutc:
            firstutc = in_spec.utc
        elif in_spec.utc > lastutc:
            lastutc = in_spec.utc
        nave = nave + 1
        ave_spec.ydataA = ave_spec.ydataA + (in_spec.ydataA * in_spec.durationSec)
        # keep track of observing time for weighted sum
        ave_spec.durationSec = ave_spec.durationSec + in_spec.durationSec

    return ave_spec, nave, firstutc, lastutc

def normalize_spec( ave_spec, firstutc, lastutc):
    """
    Normaize the average after completing sum
    input/output
    ave_spec   input raw sum of observation
               output normalized sum of observations
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
    return ave_spec

def read_hot( names, ave_hot, doScaleAve):
    """
    read_hot() reads in all files in the names list and averages hot load
    observations.   The hot load file is returned.
    While reading the files, minmum elevation and galactic latitudes are tracked
    """

    rs = radioastronomy.Spectrum()
    nhot = 0       # init count of hot files
    minGlat = +90.
    maxGlat = -90.
    minel = 200.
    maxel = -200.

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
            print(( 'File is not an astronomy file: %s' % (filename)))
            continue
        extension = parts[nparts-1]
        extension = extension.upper()
        if extension not in ['HOT', 'AST', 'CLD']:
            print(( 'Extension not recognized : %s' % ( parts[nparts-1])))
            continue

        rs.read_spec_ast(filename)
        if doScaleAve:
            rs.ydataA = rs.ydataA/rs.count
        else:
            rs.ydataA = rs.ydataA/rs.nave

        nChan = len( rs.ydataA)
        if nChan not in (32, 64, 128, 256, 512, 1024, 2048, 4096, 8192):
            print("Unusual data length (%d) for file %s" % (nChan, filename))
            print("Skipping use in averages")
            continue
        rs.azel2radec()    # compute ra,dec from az,el

        # if a hot load observation
        if rs.telel < 0:
            if nhot == 0:
                firstutc = rs.utc
                lastutc = rs.utc
        # accumulate spectra, nhot is updated also
            ave_hot, nhot, firstutc, lastutc = \
                average_spec( ave_hot, rs, nhot, firstutc, lastutc)
            if minel > rs.telel:
                minel = rs.telel
            if maxel < rs.telel:
                maxel = rs.telel
        else: # else above horizon, find min, max galactic latitudes
            if rs.gallat > maxGlat:
                maxGlat = rs.gallat
            if rs.gallat < minGlat:
                minGlat = rs.gallat
            if minel > rs.telel:
                minel = rs.telel
            if maxel < rs.telel:
                maxel = rs.telel
        # end for all files

    if nhot > 0:
        print(( "Found %3d Hot load observations" % (nhot)))
        ave_hot = normalize_spec( ave_hot, firstutc, lastutc)
    else:
        print( "No Hot load data, can not calibrate")
        sys.exit()
    return ave_hot, minel, maxel, minGlat, maxGlat
    # end of read_hot
    
def read_angles( names, lowel):
    """
    read_angles() reads all files and counts number of files with
    high elevation and high galactic latitude
    Inputs:
       lowel   minimum elevation to accept for cold load
    """

    # flag starting a new sum of cold (high elevation and galactic latitude) obs
    ncold = 0
    rs = radioastronomy.Spectrum()
    nName = len(names)
    minel = 90.
    maxel = -90.
    minGlat = 90
    maxGlat = -90
    # now average coldest data for calibration
    for filename in names:

        rs.read_spec_ast(filename)

        rs.azel2radec()    # compute ra,dec from az,el

        if rs.telel < lowel:  #if elevation too low for a cold load obs
            continue
        if rs.gallat > maxGlat:
            maxGlat = rs.gallat
        if rs.gallat < minGlat:
            minGlat = rs.gallat
        if minel > rs.telel:
            minel = rs.telel
        if maxel < rs.telel:
            maxel = rs.telel

        # count number of high elevation files
        if rs.telel > lowel:
            ncold = ncold + 1

    print("Found %d high elevation spectra in %d files" % (ncold, nName))
    return ncold, minel, maxel, minGlat, maxGlat

def velocity_indecies( in_spec, nrms, minvel, maxvel):
    """
    Compute indicies corresponding to min and max velocity
    inputs in_spec spectrum in radio astronomy structure
    nuRefFreq - reference frequency (Hz)
    nrms - number of channels to include for indicities
    minvel - minimum velocity for inner channel of lower velocity range (km/s)
    maxvel - maximum velocity for inner channel of highher velocity rane (km/s)
    """
    xv = in_spec.xdata
    nData = len( xv)        # get data length and create velocity array
    vel = np.zeros(nData)

    # create index array
    for jjj in range (0, nData):
        vel[jjj] = C * (in_spec.refFreqHz - xv[jjj])/in_spec.refFreqHz

    # xa and xb are the indices to the velocity array closest
    # to the minimum and maximum velocities
    xa, xb = gf.velocity_to_indicies( vel, minvel, maxvel)
    #try to use the channels just outside the velocity ranges defined
    nrms = 25 # define the number of channels for RMS calculation
    xa0 = xa - nrms
    if xa0 < 1:   # if before end beging of spectra, must use first channels
        xa0 = 1
    xa = xa0 + nrms
    xbe = xb + nrms
    if xbe > nData - 1:  # if after end of data must use end channels
        xbe = nData - 1
        xb = nData - nrms

    return vel, xa0, xa, xb, xbe

def read_cold( names, ave_cold, lowel, lowGlat, doScaleAve):
    """
    read_cold() reads all files, averages selected files
    (high elevation and galactic Latitude).
    Inputs:
       lowel   minimum elevation to accept for cold load
       lowGlat minimum galactic latitude
       doScaleAve - rescale counts in ascii files to match real time plotsf
    """

    # flag starting a new sum of cold (high elevation and galactic latitude) obs
    ncold = 0
    rs = radioastronomy.Spectrum()
    nName = len(names)
    # now average coldest data for calibration
    for filename in names:

        rs.read_spec_ast(filename)
        if doScaleAve:
            rs.ydataA = rs.ydataA/rs.count
        else:
            rs.ydataA = rs.ydataA/rs.nave
        rs.azel2radec()    # compute ra,dec from az,el

        if rs.telel < lowel:  #if elevation too low for a cold load obs
            continue

        # note this test excludes low galactic latitude ranges
        if rs.gallat > lowGlat or rs.gallat < -lowGlat:
            if ncold == 0:
                firstutc = rs.utc
                lastutc = rs.utc
            ave_cold, ncold, firstutc, lastutc = \
                average_spec( ave_cold, rs, ncold, firstutc, lastutc)

    # end of all files to average
    if ncold < 1:
        print( "No high elevation data: can not calibrate")
        sys.exit()
    else:
        ave_cold = normalize_spec( ave_cold, firstutc, lastutc)

    return ave_cold, ncold

def compute_gain( hv, cv, xa0, xa, xb, xbe, thot, tcold):
    """
    Compute gain on a channel by channel basis for tRx calculation
    hv, cv are arrays of hot and cold intensities in counts
    xa0, xa  are the channel ranges for one tsys measurement
    xb, xbe are the channel ranges for the other tsys measurement
    thot, tcold are hot and cold load estimates (Kelvins)

    """
    nData = len(hv)
    gainHC = np.zeros(nData)
    for iii in range(nData):
        gainHC[iii] = (hv[iii] - cv[iii])/(thot - tcold)
        if gainHC[iii] < EPSILON:
            gainHC[iii] = EPSILON

    # the hot/cold gain ratios are only used to compute tRxMiddle
    trx = np.zeros(nData)
    for iii in range(nData):
        trx[iii] = (cv[iii]/gainHC[iii]) - tcold

    tRxA = np.median(trx[xa0:xa])
    tRxB = np.median(trx[xb:xbe])
    tRxMiddle = (tRxA + tRxB)*.5

    tStdA = np.std(trx[xa0:xa])
    tStdB = np.std(trx[xb:xbe])
    tRms  = (tStdA + tStdB) * .5

    #after this point, only the hot load observations are used to compute sys
    #No cold observations are used for calibration, except for computing tRxMiddle

#    print("Median Receiver Temp: %7.3f +/- %6.3f (%6.3f %6.3f) (K)" % \
#          (tRxMiddle, tRms, tStdA, tStdB))

    # for remainder of calculations only use hot counts for calibration
    # Using hot load only reduces interference effects
    gain = np.zeros(nData)
    for iii in range(nData):
        gain[iii] = hv[iii]/(thot + tRxMiddle)
        if gain[iii] < EPSILON:
            gain[iii] = EPSILON

    gainA = np.median(gain[xa0:xa])
    gainB = np.median(gain[xb:xbe])
    gainAve = 2.0/(gainA + gainB)  # Report gain in K per Count

    return gain, gainAve, tRxMiddle, tRms, tStdA, tStdB

def flagCenter( tsky, nData):
    """
    interpolate over the central channel of the spectrum
    """

    icenter = int(nData/2)
    tsky[icenter] = (tsky[icenter-2] + tsky[icenter+2])*.5
    tsky[icenter-1] = (3.*tsky[icenter-2] + tsky[icenter+2])*.25
    tsky[icenter+1] = (tsky[icenter-2] + 3.*tsky[icenter+2])*.25
    return tsky

# code for keeping hot and cold avefrage spectra
def keep_hotcold( ave_hot, ave_cold, cpuIndex, directory):
    """
    write the hot and cold file files out
    """
    coldname = radioastronomy.utcToName( ave_cold.utc)
    coldname = coldname + ".ast"  # output in counts
    # add telescope index
    coldname = ("T%d-" % cpuIndex)+coldname
    ave_cold.ydataA = ave_cold.ydataA * ave_cold.nave
    ave_cold.write_ascii_file(directory, coldname)
    print( "Wrote Average Cold Load File: %s%s" % (directory, coldname))
    # now hot file
    hotname = radioastronomy.utcToName( ave_hot.utc)
    hotname = hotname + ".hot"  # output in counts
    # add telescope index
    hotname = ("T%d-" % cpuIndex)+hotname
    ave_hot.ydataA = ave_hot.ydataA * ave_hot.nave
    ave_hot.write_ascii_file(directory, hotname)
    print( "Wrote Average Hot  Load File: %s%s" % (directory, hotname))

    return hotname, coldname
    # end of keep_hotcold()

def tsky_gain( yv, gain):
    """
    Compute TSky based on hot and coold load observations
    Note that here the yv[], hv[] and cv[] values are normalized by the
    total number of integrations.  The raw data, read in, are increased
    linearly by the number of spectra averaged.
    """
    nData = len(yv)
    tsky  = np.zeros(nData)    # initialize arrays

    for jjj in range (0, nData):
        tsky[jjj] = yv[jjj]/gain[jjj]

    return tsky
# end of tsky_gain()
