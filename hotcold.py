#Python Scripts to read a list of files and compute raw hot and cold average spectra
#HISTORY
#18DEC10 GIL initial version from m.py
#
import sys
import datetime
import numpy as np
import radioastronomy
import copy
import statistics

EPSILON = 1.
clight = 299792458. # speed of light in m/sec

def hotaverage( names):
    """
    Read in all spectra from a list of names and return the average hot spectrum
    """
    rs = radioastronomy.Spectrum()   # create input and average structures
    hot = radioastronomy.Spectrum()
    nhot = 0

    # for all input files
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
            print 'File is not an astronomy file: ',filename
            continue
        else:
            extension = parts[nparts-1]

        extension = extension.upper()
        if extension != 'HOT':  # speed up by only  looking at hot load files
            continue
        
        rs.read_spec_ast(filename)

        if rs.telel > 0: # only working with hot load, skip elevation > 0.
            continue
            
        if nhot == 0:
            hot = copy.deepcopy( rs)
            nhot = 1
        else:
            hot.ydataA = hot.ydataA + rs.ydataA
            hot.count = hot.count + rs.count
            hot.durationSec = hot.durationSec + rs.durationSec
            nhot = nhot + 1
            
    if nhot > 0:
        hot.ydataA = hot.ydataA / float(nhot)
    else:
        print "No hot load data, can not calibrate"
            
    return nhot, hot

# end of def hotaverage

def coldaverage( names):
    """
    Compute a reference cold spectrum from a list of input file names
    """

    rs = radioastronomy.Spectrum()   # create input and average structures
    cold = radioastronomy.Spectrum()

# assume only a limited range of galactic latitudes are available
# not range above +/-60.
    use60Range = False
    minGlat = 90.                    # initialize to extremea
    maxGlat = -90.
    maxEl = -90.
    ncold = 0

    # for all input files
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
            print 'File is not an astronomy file: ',filename
            continue
        else:
            extension = parts[nparts-1]

        extension = extension.upper()
        if extension != 'AST':  # speed up by only looking at astronomy files
            continue
        
        rs.read_spec_ast(filename)

        if rs.telel < 0: # only working with observations, skip elevation <= 0.
            continue

        maxGlat = max( rs.gallat, maxGlat)
        minGlat = min( rs.gallat, minGlat)
        maxEl = min( rs.telel, maxEl)
    # finish first run through of file names; now check min and max latitudes

    # if any high galactic latitudes, use only above +/-60d 
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

    # only use the elevations above 60 degrees, if any
    if maxEl > 60.:
        maxEl = 60.
    else:
        maxEl = maxEl - 10.    #else must use highest elevations available

    # now average coldest data for calibration
    for filename in names:

        rs.read_spec_ast(filename)
        rs.azel2radec()    # compute ra,dec from az,el

        if rs.telel < maxEl:
            continue

        if rs.gallat > maxGlat or rs.gallat < minGlat:
            if ncold == 0:
                cold = copy.deepcopy( rs)
                cold.ydataA = rs.ydataA
                ncold = 1
            else:
                cold.ydataA = cold.ydataA + rs.ydataA
                cold.count = cold.count + rs.count
                cold.durationSec = cold.durationSec + rs.durationSec
                ncold = ncold + 1

    #end for all high latitude sources
    if ncold < 1.:
        print "No high galactic latitude data: can not calibrate"
    else:
        cold.ydataA = cold.ydataA/float(ncold)

    return ncold, cold
# end of def coldaverage

def compute_tsky_hot( xv, yv, hv, thot, tcold):
    """
    compute_tsky() computes the system temperature if only hot load data are available
    """

    nData = len(yv)    
    epsilons = np.full( nData, EPSILON)
    tsys = np.zeros(nData)    # initialize arrays

    Z = np.zeros(nData)
    oneMZ = np.zeros(nData)
    # For full Temp calibration, a spectrum taken at high elevation away from 
    # The galactic plan is used.   For this program the cold spectrum must be
    # the spectrum being calibrated.   See the M command for comparision
    epsilons = np.full( nData, EPSILON)
    yv = np.maximum( yv, epsilons)
    hv = np.maximum( hv, epsilons)
        # comput the cold/hot ratio
    Z = yv/hv
    oneMZ = np.full( nData, 1.) - Z
    oneMZ = np.maximum( oneMZ, epsilons)

    # the cold, receiver, temperature is this function
    tsys = ((Z*thot) - tcold)/oneMZ
    
    n6 = int(nData/6)
    n56 = 5*n6

    tsysmedian = statistics.median( tsys[n6:n56])

    tsky  = np.zeros(nData)    # initialize arrays
    S     = np.zeros(nData)    # initialize arrays

    # The system gain S is computed assuming a tsys is the cold load
    S = np.full( nData, tsysmedian+thot)/hv
    # scale the observed instensity in counts to Kelvins.
    tsky = S*yv

    return tsky
#nd of compute_tsky_hot

def compute_gain( hv, cv, thot, tcold):
    """
    compute_tsky() computes the system temperature if only hot load data are available
    """

    nData = len(cv)    
    epsilons = np.full( nData, EPSILON)

    # For full Temp calibration, a spectrum taken at high elevation, away from 
    # The galactic plan is used.
    dv = hv - cv 
    dv = np.maximum( dv, epsilons)
    hv = np.maximum( hv, epsilons)

    # use the comparison of hot and cold to get gain without trx contribution
    gain = (thot - tcold)/dv

    trx = gain*cv 
    trx = trx - tcold     # correct for cold load input in counts

    n6 = int(nData/6)
    n56 = 5*n6

    trxmedian = statistics.median( trx[n6:n56])

    print "Trx: ",trxmedian
    gain = (thot + trx)/hv
    
    return trxmedian, gain # channel by channel gain in K/counts
#end of compute_gain

def compute_tsky( xv, yv, gain, nureference):

    tsky = yv * gain

    return tsky

def compute_baseline( xv, yv, xa, xb):

    ya = statistics.median(yv[(xa-10):(xa+10)])
    yb = statistics.median(yv[(xb-10):(xb+10)])
    slope = (yb-ya)/(xb-xa)
#                baseline = tsky
#                print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
    for iii in range( nData):
        yv[iii] = (ya + (slope*(iii-xa)))
    
    return yv
# end of compute_baseline

def chan2freq( rs, chan, nureference):
    """
    Compute frequencies corresponding to the input list of channels
    """
    dx = rs.bandwidthHz/float(rs.refChan)
    dchan = float( chan - rs.refChan)
    freq = rs.centerFreqHz + (dx*dchan)
    return freq
# end of chan2freq


def flagCenter( yv):             # if flagging spike in center of plot
    """ 
    flagCenter: interpolates over the very center channels of an array
    """
    nData = len(yv)    
    icenter = int(nData/2)
    yv[icenter] = (yv[icenter-2] + yv[icenter+2])*.5
    yv[icenter-1] = (3.*yv[icenter-2] + yv[icenter+2])*.25
    yv[icenter+1] = (yv[icenter-2] + 3.*yv[icenter+2])*.25
    return yv
# end of flagCenter

