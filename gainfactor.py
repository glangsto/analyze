#Python function to compute gain factor for different processor indicies, dates and elevations
#HISTORY
#20APR06 GIL initial version based on holdcold.py
#
import sys
import datetime
import numpy as np

# define a few parameters
NMEAS = 10
dates = np.zeros(NMEAS)
els = np.zeros(NMEAS)
nMeas = np.zeros(NMEAS)
nObs = 1
nProcessors = 8
nMeasurements = 3
measurements = np.zeros((nObs, nProcessors, nMeasurements))
# now fill in measuremnts for this observations
iObs = 0
els[iObs] = 90.
nMeas[iObs] = 3
measurements[iObs, 2,] = ( 26.3, 26.0, 25.3)
measurements[iObs, 3,] = ( 49.6, 47.7, 44.9)
measurements[iObs, 4,] = ( 47.2, 48.2, 48.4)
measurements[iObs, 5,] = ( 40.0, 39.6, 38.1)
gainFactors = np.zeros((nObs, nProcessors))
normalized = False

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

def compute_gain_factor( pIndex, aveutc, el):
    """
    Compute a processor based gain factor for inputs
    Inputs:
    pIndex          Index to processor based measurements
    aveutc          Date of observation for which the factors are measured
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

