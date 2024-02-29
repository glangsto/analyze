"""
Read in an observation summary and fit the times of galaxy crossings.
From these measurements estimate the Azimuth and Elevation of the
telescope pointing.  Then compute the azimuth and elevation offsets
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 24Feb27 GIL double check the dAz offset, add cos(el) factor
# 24Feb19 GIL major reorganization fitting two crossings separately
# 24Feb19 GIL update for 2024, check for not finding offset
# 23Oct13 GIL check Az offset with FWHM galactic intensities
# 23Oct12 GIL count number of observations at each elevation
# 23Oct11 GIL fix limits to 1 day of drift scan, add log of offset
# 23Oct10 GIL limit fit to 1 days observations
# 23Apr01 GIL test using 1/2 maximum times to get dAz
# 23Mar31 GIL check dAz
# 22May06 GIL iteratively fit peaks
# 22May04 GIL get azimuth, elevation and telescope position
# 22May02 GIL initial version

import sys
import os
import datetime
import numpy as np
from matplotlib import pyplot as plt
from pylab import *
from scipy.optimize import curve_fit
from vmedian import vmedian
from vmedian import vsmooth

import gainfactor as gf
import radioastronomy

# ine file paths
#offsetfilename = "/home/karl/Research/analyze/GalacticPlaneCrossingOffsets.txt"
offsetfilename = "/Research/analyze/crossings.txt"
dazdelfilename = "/home/karl/2024-DazDel.txt"
doTwo = True
doThree = True
doFour = True

nPoleRa = (12+(51.4/60.))*15 # (degrees)
nPoleDec = (12+(51.4/60.))*15 # (degrees)
sPoleRa = nPoleRa - 180.
sPoleDec = -nPoleDec

def readGalacticOffsets( offsetfilename):
    """
    Read a table of pre-calculated Galactic lat=0 RA and Dec positions
    input - offset file name
    outputs
    decs - array of declinations for which Galactic Plane is crossed
    ra1s - array of 1st galactic plane crossing RAs
    ra2s - array of 2nd galactic plan crossing RAs
    dRas - array of deltas in RA crossings
    """
    decs = []
    ra1s = []
    ra2s = []
    dRas = []
    n = 0
    # get the home directory
    home = os.path.expanduser("~")
    f = open( home+offsetfilename, "r")
    while True:
        aline = f.readline()
        n = n + 1
        if len(aline) < 1:
            if n > 5:    # allow no more than a few blank lines
                break
            else:
                continue
        aline = aline.strip()
        if len(aline) < 1:
            continue
        if aline[0] == "#":
            continue
        parts = aline.split()
        nparts = len(parts)
        if nparts < 1:
            break
        if nparts >= 4:
            decs.append( float(parts[0]))
            ra1s.append( float(parts[1]))
            ra2s.append( float(parts[2]))
            dRas.append( float(parts[3]))
        else:
            print("Error Parsing Galactic Offsets line:")
            print(parts)

    decs = np.asarray( decs)
    ra1s = np.asarray( ra1s)
    ra2s = np.asarray( ra2s)
    n1 = len(ra1s);
    for iii in range(n1):
        if ra1s[iii] < 0.:
            ra1s[iii] = 360. + ra1s[iii]
        if ra1s[iii] > 360.:
            ra1s[iii] = ra1s[iii] - 360.
    n2 = len(ra2s)
    for iii in range(n2):
        if ra2s[iii] < 0.:
            ra2s[iii] = 360. + ra2s[iii]
        if ra2s[iii] > 360.:
            ra2s[iii] = ra2s[iii] - 360.
            
            
    dRas = np.asarray( dRas)

    f.close()

    # end of readGalacticOffsets
    return decs, ra1s, ra2s, dRas

def getCrossingsFromDec( dec, decs, ra1s, ra2s):
    """
    get the closest 1st and 2nd Ra crossings, for the input dec
    where 
    dec is deduced from observations,  
    decs, ra1s, and ra2s are read from a precomputed file.
    """

    n = len(ra1s)
    iDec = 0
    minDec = dec - decs[iDec]
    # prepare to find offsets
    if minDec < 0:
        minDec = - minDec

    # find closest dec.  Later may interpolate       
    for iii in range(n):
        dDec = dec - decs[iii]
        if dDec < 0:
            dDec = - dDec
        if dDec < minDec:
            iDec = iii
            minDec = dDec

    ra1 = ra1s[iDec]  # return pre-computed crosssing for this  declination.
    ra2 = ra2s[iDec]

    return ra1, ra2   # end of getCrossingFromDec
            
def writeDazDel( dazdelfilename, utc1, utc2, cpuIndex, az, el, dAz, dEl):
    """
    appends a new measurement of the Azimuth and elevation offsets
    inputs
    dazdelfile - file to append measuremnts to
    representative time of measurement
    az, el - recorded azimuth and elevation of observations
    daz, del - offsets to be added to az, el to get true az,el of observation
    """
    if os.path.exists( dazdelfilename):
        f = open( dazdelfilename, "a")
    else:
        f = open( dazdelfilename, "w")

    utc1str = ("%s" % utc1)        # convert utc to str
    utcparts = utc1str.split(".")  # separate by period in seconds
    utc1str = utcparts[0]          # only print seconds, no fraction
    utc2str = ("%s" % utc2)        # convert utc to str
    utcparts = utc2str.split(".")  # separate by period in seconds
    utc2str = utcparts[0]          # only print seconds, no fraction

    outline =  "%s %s %2d %8.2f %8.2f %6.2f %6.2f \r\n" % \
        (utc1str, utc2str, cpuIndex, az, el, dAz, dEl)

    f.write(outline)
    f.close()
    # end off writeDazDel()
    return

def decFromCrossingDelta( dRa, gDecs, gdRas):
    """
    decFromCrossingDelta returns the closest Declination
    matching the measured offset between two galactic crossing positions
    input
    dRa - measured right ascension difference between two galactic crossings
    gDecs - array of declinations for input offsets
    gRa1s - 1st right ascension crossing point for each declination
    dRa2s - 2nd right ascension crossing point for each declination
    gdRas - difference between ra1s and ra2s for each declination
    """
    idec = 0
    n = len( gdRas)
    # find the minimum difference in dRa
    minDRa = np.abs( gdRas[idec] - dRa)
    iminDec = idec

    print("decFrom input dRa: %7.2f" % (dRa))
    
    # find dec of crossing point separation
    for idec in range(n):
        dd = np.abs( gdRas[idec] - dRa)
        if dd < minDRa:
            iminDec = idec;
            minDRa = dd

    print( "dRa = %7.2f (%7.2f) coresponds to dec %7.1f (%d)" % \
           (dRa, gdRas[iminDec], gDecs[iminDec], iminDec))
    return gDecs[iminDec]

def gauss(x,mu,sigma,A):
    """
    Return the intensity for a single gaussian model
    input order is x axis values
    mu = center postion
    sigma = width of gaussian
    A = amplitude of gaussian
    """
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    """
    Return the intensity for a double gaussian model
    """
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def sortParams( inparams, sigmas):
    """
    sortParams re-orders parameters in order of peak height
    the params are sets of three values: time, width and peak
    returns sorted copy of orginal
    """
    n = len( inparams)
    npeaks = int(n/3)
    # for all pairs of peaks
    for i in range(npeaks-1):
        for j in range(i+1,npeaks):
            peaka = inparams[(i*3) + 2]
            peakb = inparams[(j*3) + 2]
            # if peak out of order
            if peakb > peaka:
                # swap three sets of values
                temp0 = inparams[(i*3)+0]
                temp1 = inparams[(i*3)+1]
                temp2 = inparams[(i*3)+2]
                sigma0 = sigmas[(i*3)+0]
                sigma1 = sigmas[(i*3)+1]
                sigma2 = sigmas[(i*3)+2]
                inparams[(i*3)+0] = inparams[(j*3)+0]
                inparams[(i*3)+1] = inparams[(j*3)+1]
                inparams[(i*3)+2] = inparams[(j*3)+2]
                inparams[(j*3)+0]  = temp0
                inparams[(j*3)+1]  = temp1
                inparams[(j*3)+2]  = temp2
                sigmas[(j*3)+0]  = sigma0
                sigmas[(j*3)+1]  = sigma1
                sigmas[(j*3)+2]  = sigma2
    return inparams, sigmas

# for diagnostics
doTest = False

def selectIntensities( filename, gDecs, gRa1s, gRa2s, gdRas):
    """
    selectIntensities() reads a "T" integrated intensity file and
    selects only observations for the elevation with the most measurements
    The function returns the selected intensities
    """

    # read date (string), offset seconds
    firstdate, utcIns, secIns, tSumIns, dTIns, azs, els, inRas, inDecs, \
        gallons, gallats = gf.readAllValues( filename)
    nData = len(secIns)
        
    ellist = np.zeros(100)   # assume up to 100 different elevations
    elcount = np.zeros(100)
    elDec = np.zeros(100)
    azlist = np.zeros(100)
    nel = 0
    newEl = True
    
    for iii in range(nData):
        el = els[iii]
        for ii in range(nel):
            if ellist[ii] == el:
                elcount[ii] = elcount[ii] + 1
                newEl = False
        if newEl:         # if el not in list
            ellist[nel] = el  #new el
            elcount[nel] = 1
            elDec[nel] = inDecs[iii]
            azlist[nel] = azs[iii]
            nels = nel + 1

    #report number of els
    print("Found %d elevations %f ... %f" % \
          (nels, ellist[0], ellist[nels-1]))
    for ii in range(nels):
        print("%d: el: %f count: %d" % \
          (ii+1, ellist[ii], elcount[ii]))

    maxElCount = 0
    maxElI = 0
    # now find the elevation with most measurements.
    for ii in range(nels):
        if maxElCount < elcount[ii]:
            maxElI = ii
            maxElCount = elcount[ii]
            maxDec = elDec[ii] 
    print("Max measurements elevation: %.1f (d) = Dec: %.1f" % \
          (ellist[maxElI], maxDec))
    azKeep = azlist[maxElI]
    print("Azimuth Max measurements  : %.1f (d)" % \
          (azKeep))
    decKeep = maxDec
    
    # now use first declination to select the RAs of 1st and 2nd crossings
    # these values are from the input file
    raX1, raX2 = getCrossingsFromDec( maxDec, gDecs, gRa1s, gRa2s)    

    print("Ra crossing points: %.1f, %.1f" %  (raX1, raX2))
    # also assume els and ras can not be off by more than a few degrees
    fitRaRange = 12.  # assume must be close to actual
    raX1min = raX1 - fitRaRange
    raX1max = raX1 + fitRaRange
    raX2min = raX2 - fitRaRange
    raX2max = raX2 + fitRaRange

    print("Ra range 1: %.1f, %.1f" %  (raX1min, raX1max))
    print("Ra range 2: %.1f, %.1f" %  (raX2min, raX2max))

    # remember the az,el with the obs with max number of elevation samples
    elKeep = ellist[maxElI]

    # prepare to save different ra ranges
    inRa1s = np.zeros(nData)
    tSum1s = np.zeros(nData)
    inRa2s = np.zeros(nData)
    tSum2s = np.zeros(nData)

    # divid the data into two parts separated by galactic north pole
    nRa1 = 0
    nRa2 = 0

    for ii in range(nData):
        if els[ii] == elKeep:
            ra = inRas[ii]
            # if found an Ra in the 1st range
            if not (ra < nPoleRa and ra > sPoleRa):
                continue
            inRa1s[nRa1] = ra
            tSum1s[nRa1] = tSumIns[ii]
            nRa1 = nRa1 + 1

    for ii in range(nData):
        if els[ii] == elKeep:
            ra = inRas[ii]
            # if found an Ra in the 2nd range
            if ra < nPoleRa and ra > sPoleRa:
                continue
            inRa2s[nRa2] = ra
            tSum2s[nRa2] = tSumIns[ii]
            nRa2 = nRa2 + 1

    # trim to only input values
    inRa1s = inRa1s[0:nRa1]
    tSum1s = tSum1s[0:nRa1]
    inRa2s = inRa2s[0:nRa2]
    tSum2s = tSum2s[0:nRa2]

    iiRa1 = 0
    iiRa2 = 0
    
    # now further trim data to those to be fit
    for ii in range(nRa1):
        if inRa1s[ii] > raX1min and inRa1s[ii] < raX1max:
            inRa1s[iiRa1] = inRa1s[ii]
            tSum1s[iiRa1] = tSum1s[ii]
            iiRa1 = iiRa1 + 1
        
    # now further trim data to those to be fit
    for ii in range(nRa2):
        if inRa2s[ii] > raX2min and inRa2s[ii] < raX2max:
            inRa2s[iiRa2] = inRa2s[ii]
            tSum2s[iiRa2] = tSum2s[ii]
            iiRa2 = iiRa2 + 1

    print("Selected %d out of %d messurements for 1st crossing" % \
          (iiRa1, nRa1))
    print("Selected %d out of %d messurements for 2nd crossing" % \
          (iiRa2, nRa2))

    nRa1 = iiRa1
    nRa2 = iiRa2

    # trim again to use only values near the peaks
    inRa1s = inRa1s[0:nRa1]
    tSum1s = tSum1s[0:nRa1]
    inRa2s = inRa2s[0:nRa2]
    tSum2s = tSum2s[0:nRa2]

    if doTest:
        print("Found %3d samples in Ra range %5.1f to %5.1f" % \
              (nRa1, raX1min, raX1max))
        print("Found %3d samples in Ra range %5.1f to %5.1f" % \
              (nRa2, raX2min, raX2max))

    # now find maximium intensity in eacsh range
    maxi1 = 0
    maxi2 = 0
    for iii in range(nRa1):
        if tSum1s[maxi1] < tSum1s[iii]:
            maxi1 = iii
    for iii in range(nRa2):
        if tSum2s[maxi2] < tSum2s[iii]:
            maxi2 = iii
        
    if doTest:
        print("Found Crossing max 1 at %5.1f, %6.1f" % \
              (inRa1s[maxi1], tSum1s[maxi1]))
        print("Found Crossing max 2 at %5.1f, %6.1f" % \
              (inRa2s[maxi2], tSum2s[maxi2]))

    if doTest:
        plt.plot(inRas, tSumIns, color='blue',lw=3,
                 label='Intensities')
        plt.plot(inRa1s, tSum1s, color='cyan',lw=3,
                 label='Crossing 1')
        plt.plot(inRa2s, tSum2s, color='red',lw=3,
                 label='Crossing 2')
        plt.show()

    # end of selectIntensities()
    return firstdate, azKeep, elKeep, decKeep, nRa1, inRa1s, tSum1s, nRa2, inRa2s, tSum2s

def fitCrossing( filename, gDecs, gRa1s, gRa2s, gdRas):
    # the galacitic RAs crossings are pre-computed and written to a file
    gDecs, gRa1s, gRa2s, gdRas = readGalacticOffsets( offsetfilename)

    """
    fitCrossing takes the integrated intensities from a sumary file,
    then fits a series of gaussians
    This version fits one gaussian at a time, until the fit fails,
    then simultaneously fits all gaussians
    """

    # read in all values from T, separated into two arrays
    firstdate, azKeep, elKeep, decKeep, nRa1, ra1s, tSum1s, nRa2, ra2s, \
        tSum2s = selectIntensities( filename, gDecs, gRa1s, gRa2s, gdRas)

    nData = nRa1
    # create a temporary array to interatively fit.
    tTemps = np.zeros(nData)
    t1Max = 0.
    i1Max = 0
    for i in range(nData):
        tTemps[i] = tSum1s[i]
        if t1Max < tTemps[i]:
            t1Max = tTemps[i]
            i1Max = i
    ra1Max = ra1s[i1Max]
    
    # set number of gaussians to fit and init array of peak values
    iMaxs = [i1Max, 0, 0, 0, 0]
    tMaxs = [t1Max, 0., 0., 0., 0.]
    # keep fit results
    ra1Peaks = [ra1Max, ra1Max, 0., 0., 0.]
    dRa1Peaks = [0., 0., 0., 0., 0.]
    tPeaks = [t1Max, t1Max, 0., 0., 0.]
    dTPeaks = [0., 0., 0., 0., 0.]
    widths = [50., 50., 0., 0., 0.]
    dWidths = [0., 0., 0., 0., 0.]

    # estimate the default width of gaussian in RA direction (degrees)
    width = 50.
    # Two crossings can be measured only if below minimum declination 
    maxDecCrossing = 55.

    ra1Min = 360.
    ra1Max = 0.
    # now find the range of data for first crossing
    for iii in range(nData):
        if ra1Min > ra1s[iii]:
            ra1Min = ra1s[iii]
        if ra1Max < ra1s[iii]:
            ra1Max = ra1s[iii]

    dT = 20.

    # get the indexs to  the first intensity peak
    i1Max = np.argmax( tTemps)
    t1Max = tTemps[i1Max]
    r1Max = ra1s[i1Max]
    ng = 0

    # estimate the fit the gaussian
    expected=(r1Max, width, t1Max)

    bestFit = ng
    sigma1 = [ 10., 10., dT]
    try:
#        params1,cov1=curve_fit(gauss,ra1s,tTemps,expected,sigma=sigma1)
        params1,cov1=curve_fit(gauss,ra1s,tTemps,expected)
        sigma1=sqrt(diag(cov1))
        bestFit = ng
    except:
        print("Error fitting gaussian %d" % (ng+1))
        params1 = expected

    print("Fit %2d:  %9.2f %8.2f %7.2f" % \
          (1, params1[0], params1[2], params1[1]))
    print(" +/-  :  %9.2f %8.2f %7.2f" % \
          (sigma1[0], sigma1[2], sigma1[1]))

# now start fitting 2nd crossing
    nData = nRa2
    ra2Min = 360.
    ra2Max = 0.
    for iii in range(nData):
        if ra2Min > ra2s[iii]:
            ra2Min = ra2s[iii]
        if ra2Max < ra2s[iii]:
            ra2Max = ra2s[iii]
            
    # now find the range of data for first crossing
    # create a temporary array to interatively fit.
    tTemps = np.zeros(nData)
    t2Max = 0
    i2Max = 0
    for i in range(nData):
        tTemps[i] = tSum2s[i]
        if t2Max < tTemps[i]:
            t2Max = tTemps[i]
            i2Max = i
    
    # get the indexs to  the first intensity peak
    i2Max = np.argmax( tTemps)
    t2Max = tTemps[i2Max]
    ra2Max = ra2s[i2Max]

    # estimate the fit the gaussian
    expected=(ra2Max, width, t2Max)

    bestFit = ng
    sigma2 = [ 10., 10., dT]
    try:
#        params1,cov1=curve_fit(gauss,ra1s,tTemps,expected,sigma=sigma1)
        params2,cov2=curve_fit(gauss,ra2s,tTemps,expected)
        sigma2=sqrt(diag(cov1))
        bestFit = ng
    except:
        print("Error fitting gaussian %d" % (ng+1))
        params2 = expected

    print("Fit %2d:  %9.2f %8.2f %7.2f" % \
          (2, params2[0], params2[2], params2[1]))
    print(" +/-  :  %9.2f %8.2f %7.2f" % \
          (sigma2[0], sigma2[2], sigma2[1]))

    # end of fitCrossing
    return firstdate, bestFit, azKeep, elKeep, decKeep, ra1s, tSum1s, ra2s, tSum2s, params1, sigma1, params2, sigma2

def main():
    """
    Main executable for gridding astronomical data
    """

    nargs = len(sys.argv)

    if nargs < 2:
        print('Fit Crossings: find the time of crossings of peak intensities')
        print('fitCrossing <Spectrum Summary-files>')
        sys.exit()
    print( "reading %d files" % (nargs - 1))

# first read through all data and find hot load
    names = sys.argv[1:]
    names = sorted(names)

#    print(names)
    firstdate = ""

    timefmt = "%Y-%m-%d %H:%M:%S"
    rs = radioastronomy.Spectrum()
 #   print( "Reading %d files" % (len(names)))

    # Read pre-computed 2 Ra and Dec crossings from input file
    gDecs, gRa1s, gRa2s, gdRas = readGalacticOffsets( offsetfilename)

    for filename in names:
        print ("File: %s" % (filename))

        # now read the file containing intensities versus elevation and ras
        firstdate, ng, azKeep, elKeep, decKeep, ra1s, tSum1s, ra2s, tSum2s, \
            params1, sigma1, params2, sigma2 = \
                fitCrossing( filename, gDecs, gRa1s, gRa2s, gdRas)

        nData = len(ra1s)
#        print("Fit successful for %d gausians" % (ng))
        ra1 = params1[0]
        ra2 = params2[0]
        # now need to keek ra1 in lower range of ras
        if ra1 > ra2:
            temp = ra1
            ra1 = ra2
            ra2 = temp

        avera = 0.5*(ra1+ra2)
        dRa = ra2 - ra1
        dAz = 0.
        # now that we have dRa, get the actual deduced 
        foundDec = decFromCrossingDelta( dRa, gDecs, gdRas)
        dEl = decKeep - foundDec
        if decKeep > 55. or decKeep < -55.:
            dEl = 0.

        #Assume the Azimuth offset is due to the horn tilted east or west.
        #pointed south
        if azKeep > 90. and azKeep < 270.:
            dAz = avera - nPoleRa
        else:
            dAz = nPoleRa - avera
        # fix T sign cocnvention
        dAz = - dAz
        dEl = - dEl
        
        print("Ave RA: %7.3fd (%7.3fh)" %
                  (avera, avera/15.))
        print("Pol RA: %7.3fd (%7.3fh)" %
                  (nPoleRa, nPoleRa/15.))
        print("dAz: %7.3fd, dEl: %7.3fd" %
              (dAz, dEl))
        print("With T comamnd use argument:")
        print("     -O %.1f %.1f" % (dAz, dEl))

        plt.plot(ra1s, tSum1s, color='blue',lw=3,
                 label='1st Data')
        plt.plot(ra2s, tSum2s, color='blue',lw=3,
                 label='2nd Data')
        if params1[0] != 0.:
            plt.plot(ra1s, gauss(ra1s,*params1),color='cyan',lw=3,
                 label='1st Crossing')
        if params2[0] != 0.:
            plt.plot(ra2s, gauss(ra2s,*params2),color='red',lw=3,
                 label='2nd Crossing')
#        print("1 Gaussian Fit:")
        i = 0
        print("      RA     +/-   Intensity    +/-     Width   +/- ")
        print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % (
            params1[0+i*3], sigma1[0+i*3],
            params1[2+i*3], sigma1[2+i*3],
            params1[1+i*3], sigma1[1+i*3]))
        print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % (
            params2[0+i*3], sigma2[0+i*3],
            params2[2+i*3], sigma2[2+i*3],
            params2[1+i*3], sigma2[1+i*3]))
        lastId = 0
        
        plt.xlabel( "Right Ascention (Degrees) %s  dAz:%5.1fd dEl:%5.1fd" % (firstdate, dAz, dEl))
        plt.ylabel( "Peak Intensity (K)")
        plt.title("%s Galactic Integrated Intensities - Tel:%2d" % (
            firstdate, lastId))
        plt.legend()

        plt.show()
        print( "Using %d values from file %s" % ( nData, filename))

    return


if __name__ == "__main__":
    main()
