"""
Read in an observation summary and fit the times of galaxy crossings.
From these measurements estimate the Azimuth and Elevation of the
telescope pointing.  Then compute the azimuth and elevation offsets
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 23Oct12 GIL count number of each elevation
# 23Oct11 GIL fix limits to 1 day of drift scan, add log of offset
# 23Oct10 GIL limit fit to 1 days observations
# 23Apr01 GIL test using 1/2 maximum times to get dAz
# 23Mar31 GIL check dAz
# 22May06 GIL iteratively fit peaks
# 22May04 GIL get azimuth, elevation and telescope position
# 22May02 GIL initial version

import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from pylab import *
from scipy.optimize import curve_fit
import datetime
import gainfactor as gf
from vmedian import vmedian
from vmedian import vsmooth

import radioastronomy

# ine file paths
offsetfilename = "/home/karl/Research/analyze/GalacticPlaneCrossingOffsets.txt"
dazdelfilename = "/home/karl/2023-DazDel.txt"
doTwo = True
doThree = True
doFour = True

GalacticPolRa = (12+(51.4/60.))*15 # (degrees)
GalacticPolDec = 27.13            # (degrees)

def readGalacticOffsets( offsetfilename):
    """
    Read a table of RA and Dec positions
    input - offset file name
    outputs
    decs - array of declinations for which Galactic Plane is crossed
    ra1s - array of 1st  galactic plane crossing RAs
    ra2s - array of 2nd galactic plan crossing RAs
    dRas - array of deltas in RA crossings
n    """
    decs = []
    ra1s = []
    ra2s = []
    dRas = []
    count = 0
    f = open( offsetfilename, "r")
    while True:
        aline = f.readline()
        aline = aline.strip()
        if len(aline) < 1:
            count = count + 1
            if count < 10:
                continue
            else:
                break
        if aline[0] == "#":
            continue
        parts = aline.split()
        nparts = len(parts)
        if nparts < 1:
            break
        if nparts == 4:
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
    dRas = np.asarray( dRas)
    # end of readGalacticOffsets
    return decs, ra1s, ra2s, dRas

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

def decFromCrossingDelta( dRa, decs, ra1s, ra2s, dRas):
    """
    decFromCrossingDelta returns the closest Declination
    matching the measured offset between two galactic crossing positions
    input
    dRa - measured right ascension difference between two galactic crossings
    decs - array of declinations
    ra1s - 1st right ascension crossing point for each declination
    ra2s - 2nd right ascension crossing point for each declination
    dRas - difference between ra1s and ra2s for each declination
    """
    idec = 0
    n = len( dRas)
    ddRa = 360 - dRa
    while ddRa < dRas[idec] and idec < n:
        idec = idec + 1

    print( "dRa = %7.1f coresponds to dec %7.1fd (%d)" %
           (dRa, decs[idec], idec))
    return decs[idec]

def gauss(x,mu,sigma,A):
    """
    Return the intensity for a single gaussian model
    """
    return A*exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    """
    Return the intensity for a double gaussian model
    """
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def trimodal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3):
    """
    Return the intensity for a triple gaussian model
    """
    return gauss(x, mu1, sigma1, A1) + \
        gauss(x, mu2, sigma2, A2) + \
        gauss(x, mu3, sigma3, A3)

def quadmodal(x,mu1,sigma1,A1,mu2,sigma2,A2, \
              mu3,sigma3,A3, mu4,sigma4,A4):
    """
    Return the intensity for a four gaussian model
    """
    return gauss(x, mu1, sigma1, A1) + \
        gauss(x, mu2, sigma2, A2) + \
        gauss(x, mu3, sigma3, A3) + \
        gauss(x, mu4, sigma4, A4)

def fivemodal(x,mu1,sigma1,A1,mu2,sigma2,A2, \
              mu3,sigma3,A3, mu4,sigma4,A4, mu5,sigma5,A5):
    """
    Return the intensity for a five gaussian model
    """
    return gauss(x, mu1, sigma1, A1) + \
        gauss(x, mu2, sigma2, A2) + \
        gauss(x, mu3, sigma3, A3) + \
        gauss(x, mu4, sigma4, A4) + \
        gauss(x, mu5, sigma5, A5)

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

def fitCrossing( filename):
    """
    fitCrossing takes the integrated intensities from a sumary file,
    then fits a series of gaussians
    This version fits one gaussian at a time, until the fit fails,
    then simultaneously fits all gaussians
    """

    # read date (string), offset seconds
    firstdate, utcIns, secIns, tSumIns, dTIns, azs, els, gallonIns, gallatIns = \
        gf.readAllValues( filename)
    # count the number of different elevations
    nels = 1
    ellist = np.zeros(100)
    elcount = np.zeros(100)
    ellist[0] = els[0]     # start with first el
    for el in els:
        newEl = True      # assume a new el
        for ii in range(nels):
            if ellist[ii] == el:
                elcount[ii] = elcount[ii] + 1
                newEl = False
        if newEl:         # if el not in list
            ellist[nels] = el  #new el
            elcount[nels] = 1
            nels = nels + 1

    #report number of els
    print("Found %d elevations %f ... %f" % \
          (nels, ellist[0], ellist[nels-1]))
    for ii in range(nels):
        print("%d: el: %f count: %d" % \
          (ii+1, ellist[ii], elcount[ii]))

    maxElCount = 0
    maxElI = 0
    for ii in range(nels):
        if maxElCount < elcount[ii]:
            maxElI = ii
            maxElCount = elcount[ii]
    print("Max measurements elevation: %.1f (d)" % \
          (ellist[maxElI]))
    firstdate, utcs, secs, tSums, dTs, azs, els, gallons, gallats = \
        gf.readAllValues( filename)
    # now select only the obs with max elevations
    nData = len(secs)
    elKeep = ellist[maxElI]
    n = 0                                     # count obs kept
    for ii in range(nData):
        if els[ii] == elKeep:
            utcs[n] = utcs[ii]
            secs[n] = secs[ii]
            secIns[n] = secIns[ii]
            tSums[n] = tSums[ii]
            tSumIns[n] = tSumIns[ii]
            dTs[n] = dTs[ii]
            dTIns[n] = dTIns[ii]
            gallons[n] = gallons[ii]
            gallats[n] = gallats[ii]
            els[n] = els[ii]
            azs[n] = azs[ii]
            n = n + 1

    # now only elevatiosn to be fit are kept
    nData = n
    secs = secs[0:n]    # trim out extra values
    tSums = tSums[0:n]
    dTs = dTs[0:n]
    utc1 = utcs[0]                             # keep track of utcs for log
    utc2 = utcs[nData-1]                       #
    n2 = int(nData/2)
    nFit = nData - 1                           # assume all data are used
#    print("%s %lf, %lf, %lf " % (firstdate, secs[n2], tSums[n2], dTs[n2]))
    #next phase is to start fit out of galactic plane,  re arrange samples
    inPlane = True                             # Assume in galactic plane 
    firstN = 0
    minI = 0
    minT = tSums[minI]                            # initialize minimum
    for iii in range(nData):
        # need to fit starting out of galactic plane
        if inPlane:
            if gallats[iii] < -45. or gallats[iii] > 45.:
                firstN = iii
                # no longer looking for start of out of plane data
                inPlane = False
                print("First point out of galactic plane: %d at %7.1f d" %
                      (firstN, gallats[firstN]))
                minT = tSums[firstN]
                minI = firstN
        else:
            if secs[iii] > secs[firstN]+86400. : # if more than one day
                nFit = iii
                break;
        # gaussian fit is better with zero baseline, find minimum
            if tSums[iii] < minT:
                minT = tSums[iii]
                minI = iii

    # recompute miminum Temp from median of 10 values around min
    minI0 = minI - 10
    minI1 = minI + 10
    if minI0 < 0:
        minI = 0
    if minI1 > nData - 1:
        maxI1 = nData - 1;
    minT = np.median(tSums[minI0:minI1])

    # trim data to 1 day
    nMove = nFit - firstN    # these samples will be moved to the beginning
    nMore = firstN           # these may be moved to after first
    if nFit == nData:
        tSums = tSumIns - minT
    else:
        nMove1 = nMove + 1
        if firstN + nMove1 > nData:
            nMove1 = nMove
        secs[0:nMove1] = secIns[firstN:firstN+nMove1]
        #select temperature - minimum
        tSums[0:nMove1] = tSumIns[firstN:firstN+nMove1] - minT
        dTs[0:nMove1] = dTIns[firstN:firstN+nMove1]
        nData = nMove1
        
    # if not a full day of observations, append the skipped data to end
    if secs[nMove] < secs[0] + 86380.:
        # move assuming a day wrap
        secs[nMove:nMove+nMore] = secIns[0:nMore] + 86400.
        tSums[nMove:nMove+nMore] = tSumIns[0:nMore]
        dTs[nMove:nMove+nMore] = dTs[0:nMore]
        nData = nMove + nMore

    secs = secs[0:nData]
    tSums = tSums[0:nData]
    tSums = vmedian(tSums,1)                   # smooth out peaks
    tSums = vsmooth(tSums,1)                   # smooth out peaks
    dTs  = dTs[0:nData]/4.                     # smoothing reduces RMS
    
    t0 = secs[0]
    t1 = t0 + 86400.
    iLimit = nData - 1
    # to fit a few gaussians, must limit to exactly 1 day of measurements
    for iii in range(nData-1):
        if secs[iii] > secs[iii+1]:
            print("Time Order Error at %d: %.2f > %.2f" % \
                  (iii, secs[iii],secs[iii+1]))
            print("nFit %d, nMove %d, nMore %d, nData %d" % \
                  (nFit, nMove, nMore, nData))
        if secs[iii] > t1:
            iLimit = iii
            break
    
    # create an array of seconds
    nt = len(tSums)

    # create a temporary array to interatively fit.
    tTemps = np.zeros(nt)
    for i in range(nt):
        tTemps[i] = tSums[i]

    # set number of gaussians to fit and init array of peak values
    NGAUSS = 5
    iMaxs = [0, 0, 0, 0, 0]
    tMaxs = [0., 0., 0., 0., 0.]
    # keep fit results
    utcPeaks = [0., 0., 0., 0., 0.]
    dUtcPeaks = [0., 0., 0., 0., 0.]
    tPeaks = [0., 0., 0., 0., 0.]
    dTPeaks = [0., 0., 0., 0., 0.]
    widths = [0., 0., 0., 0., 0.]
    dWidths = [0., 0., 0., 0., 0.]
    
    # now fit all gaussians
    for ng in range(NGAUSS):
        
        # get the indexs to  the first intensity peak
        iMax = np.argmax( tTemps)
        tMax = tTemps[iMax]
        iMaxs[ng] = iMax
        tMaxs[ng] = tMax

        # limit search for half max to range near peak
        # but not beyond end of data arrays
        if iMax < nt/2:
            nseek = int(iMax*.4)
        else:
            nseek = int((nt - iMax)*.4)
        # assume width is no more than one hour
        width = 3600.
        tHalf = tMax/2.
        # now find the half width:
        for i in range(nseek):
            j = iMax - i
            k = iMax + i
            # if found the half width
            if tTemps[j] < tHalf:
                width = secs[iMax] - secs[j]
                break
            elif tTemps[k] < tHalf:
                width = secs[k] - secs[iMax]
                break
        # end searching for width
                
#        print("t%d = %8.1f, y = %8.3f, width = %7.2f " % \
#          (ng+1, secs[iMax], tMax, width))
                
        # estimate the fit the gaussian
        expected=(secs[iMax], width, tMax)

        sigma1 = [ 0., 0., 0.]
        try:
            params1,cov1=curve_fit(gauss,secs,tTemps,expected,sigma=dTs)
            sigma1=sqrt(diag(cov1))
        except:
            print("Error fitting gaussian %d" % (ng+1))
            params1 = expected

        # fit was successful, subtract fit and try again
        for i in range(nt):
            tTemps[i] = tTemps[i] - gauss( secs[i], *params1)

        utcPeaks[ng] = params1[0]
        widths[ng] = params1[1]
        tPeaks[ng] = params1[2]
        dUtcPeaks[ng] = sigma1[0]
        dWidths[ng] = sigma1[1]
        dTPeaks[ng] = sigma1[2]

        print("Fit %2d:  %9.2f %8.2f %7.2f" % \
              (ng, utcPeaks[ng], tPeaks[ng], widths[ng]))
        print(" +/-  :  %9.2f %8.2f %7.2f" % \
              (dUtcPeaks[ng], dTPeaks[ng], dWidths[ng]))
                
    # end for all gaussians
    params1 = [utcPeaks[0], widths[0], tPeaks[0]]

    # prepare to bound fits in time
    tmin = secs[0] - 10.
    tmax = secs[nt-1] + 10.

    # now try 1, 2, 3 and 4 gaussians
    # keep the largest number that fits

    expected1 = ( utcPeaks[0], widths[0], tPeaks[0])
    sigma1 = [0., 0., 0.]
    try:
        params1,cov1=curve_fit(gauss,secs,tSums,expected1,sigma=dTs)
        sigma1=sqrt(diag(cov1))
        ng = 1
    except:
        print("Error trying a 1 gaussian fit")
        params1 = [utcPeaks[0], widths[0], tPeaks[0]]

    expected2 = ( utcPeaks[0], widths[0], tPeaks[0], 
                  utcPeaks[1], widths[1], tPeaks[1])
    bounds2 =    [ (tmin, 50., 50., tmin, 50., 50.),
                   (tmax, 90000., 300000., tmax, 90000., 300000.)]
    sigma2 = [0., 0., 0., 0., 0., 0.]
    try:
        params2,cov2=curve_fit(bimodal,secs,tSums,expected2,sigma=dTs,
                               bounds=bounds2)
        sigma2=sqrt(diag(cov2))
        ng = 2
    except:
        print("Error trying a 2 gaussian fit")
        params2 = [ utcPeaks[0], widths[0], tPeaks[0], 
                  utcPeaks[1], widths[1], tPeaks[1]]

    expected3 = ( utcPeaks[0], widths[0], tPeaks[0], 
                  utcPeaks[1], widths[1], tPeaks[1], 
                  utcPeaks[1], widths[2], tPeaks[2])
    bounds3 =    [ (tmin, 50., 50., tmin, 50., 50.,
                    tmin, 50., 50.),
                   (tmax, 90000., 300000., tmax, 90000., 300000.,
                    tmax, 90000., 300000.)]
    sigma3 = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
    try:
        
        params3,cov3=curve_fit(trimodal,secs,tSums,expected3,sigma=dTs,
                               bounds=bounds3)
        sigma3=sqrt(diag(cov3))
        ng = 3
    except:
        print("Error trying a 3 gaussian fit")
        params3 = [ utcPeaks[0], widths[0], tPeaks[0], 
                  utcPeaks[1], widths[1], tPeaks[1], 
                  utcPeaks[1], widths[2], tPeaks[2]]

            
    expected4 = ( utcPeaks[0], widths[0], tPeaks[0], 
                  utcPeaks[1], widths[1], tPeaks[1], 
                  utcPeaks[2], widths[2], tPeaks[2], 
                  utcPeaks[3], widths[3], tPeaks[3])
    bounds4 =    [ (tmin, 50., 50., tmin, 50., 50.,
                    tmin, 50., 50.,  tmin, 50., 50.),
                   (tmax, 90000., 300000., tmax, 90000., 300000.,
                   tmax, 90000., 300000., tmax, 90000., 300000.)]
    sigma4 = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
    try:
        params4,cov4=curve_fit(quadmodal,secs,tSums,expected4,sigma=dTs,
                               bounds=bounds4)
        sigma4=sqrt(diag(cov4))
        ng = 4
    except:
        print("Error trying a 4 gaussian fit")
        params4 = [utcPeaks[0], widths[0], tPeaks[0], 
                  utcPeaks[1], widths[1], tPeaks[1], 
                  utcPeaks[2], widths[2], tPeaks[2], 
                  utcPeaks[3], widths[3], tPeaks[3]]

    expected5 = ( utcPeaks[0], widths[0], tPeaks[0], 
                  utcPeaks[1], widths[1], tPeaks[1], 
                  utcPeaks[2], widths[2], tPeaks[2], 
                  utcPeaks[3], widths[3], tPeaks[3],
                  utcPeaks[4], widths[4], tPeaks[4])
    bounds5 =    [ (tmin, 50., 50., tmin, 50., 50.,
                  tmin, 50., 50.,  tmin, 50., 50.,
                  tmin, 50., 50.),
                  (tmax, 90000., 300000., tmax, 90000., 300000.,
                   tmax, 90000., 300000., tmax, 90000., 300000.,
                   tmax, 90000., 300000.)]
    sigma5 = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
    try:
        params5,cov5=curve_fit(fivemodal,secs,tSums,expected5,sigma=dTs,
                               bounds=bounds5)
        sigma5=sqrt(diag(cov5))
        ng = 5
    except:
        print("Error trying a 5 gaussian fit")
        params5 = [utcPeaks[0], widths[0], tPeaks[0], 
                   utcPeaks[1], widths[1], tPeaks[1], 
                   utcPeaks[2], widths[2], tPeaks[2], 
                   utcPeaks[3], widths[3], tPeaks[3],
                   utcPeaks[4], widths[4], tPeaks[4]]

    # now sort each result in intensity order
    # no need to sort if only 2 peaks
#    params1, sigma1 = sortParams( params1, sigma1)
    params2, sigma2 = sortParams( params2, sigma2)
    params3, sigma3 = sortParams( params3, sigma3)
    params4, sigma4 = sortParams( params4, sigma4)
    params5, sigma5 = sortParams( params5, sigma5)
    
    return firstdate, utc1, utc2, azs[0], elKeep, secs, tSums, dTs, ng, params1, sigma1, params2, sigma2, params3, sigma3, params4, sigma4, params5, sigma5

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

    print(names)
    firsttime = ""
    lasttime = ""
    firstdate = ""
    lastdate = ""
    count = 0
    utcs = []

    timefmt = "%Y-%m-%d %H:%M:%S"
    rs = radioastronomy.Spectrum()
 #   print( "Reading %d files" % (len(names)))
           
    for filename in names:
        print ("File: %s" % (filename))

        # return the data and three fits
        firstdate, utca, utcb, az, el, secs, tSums, dTs, ng, params1, sigma1, params2, sigma2, params3, sigma3, params4, sigma4, params5, sigma5 = fitCrossing( filename)

        nData = len(secs)
        print("Fit successful for %d gausians" % (ng))
        lastId, lastaz, lastel = gf.lastazel()
        lastaz = az
        lastel = el
        rs.telaz = az
        rs.telel = el
        # retrieve the telescope coordinates
        tellon, tellat, telelev = gf.lonlatelev()
        rs.tellon = tellon
        rs.tellat = tellat
        rs.telelev = telelev
        print("Telescope lon, lat, elev: %7.2f %7.2f %7.1f" % 
              (rs.tellon, rs.tellat, rs.telelev))
        print("Telescope  Az,  El, id  : %7.2f %7.2f %2d" %
              (rs.telaz, rs.telel, lastId))
        utcmidnight = datetime.datetime.strptime("20" + firstdate + " 00:00:00",
                                                 timefmt)
        utc1 = utcmidnight   # create variables, for update below
        utc2 = utc1
        if ng == 5:
            utc1 = utcmidnight + datetime.timedelta( seconds = params5[0])
            utc2 = utcmidnight + datetime.timedelta( seconds = params5[3])
        elif ng == 4:
            utc1 = utcmidnight + datetime.timedelta( seconds = params4[0])
            utc2 = utcmidnight + datetime.timedelta( seconds = params4[3])
        elif params3[0] != 0.:
            utc1 = utcmidnight + datetime.timedelta( seconds = params3[0])
            utc2 = utcmidnight + datetime.timedelta( seconds = params3[3])
        elif params2[0] != 0.:
            utc1 = utcmidnight + datetime.timedelta( seconds = params2[0])
            utc2 = utcmidnight + datetime.timedelta( seconds = params2[3])
        rs.utc = utc1            

        print("Time of first  crossing: %s" % (utc1))

        rs.azel2radec()
        ra1 = rs.ra
        dec1 = rs.dec
        print("RA, Dec     of crossing: %7.3fd %7.3fd (%7.3fh)" % 
              (rs.ra, rs.dec, rs.ra/15.))

        rs.utc = utc2
        print("Time of second crossing: %s" % (utc2))
        rs.azel2radec()
        ra2 = rs.ra
        dec2 = rs.dec
        print("RA, Dec     of crossing: %7.3fd %7.3fd (%7.3fh)" % 
              (rs.ra, rs.dec, rs.ra/15.))

        dRa = ra1 - ra2
        if dRa < 0:
            dRa = - dRa
        # read in offsets vs dec
        decs, ra1s, ra2s, dRas = readGalacticOffsets( offsetfilename)
        founddec = decFromCrossingDelta( dRa, decs, ra1s, ra2s, dRas)
        # can't determine el offset if dec > 55.
        avera = (ra1 + ra2)/2.
        dAz = 0
        if rs.dec > 55.:
            dEl = 0.
            if ng == 5:
                if params5[2] * 2. * params5[5]:
                    avera = params5[2]
            elif ng == 4:
                if params4[2] * 2. * params4[5]:
                    avera = params4[2]
            elif ng == 3:
                if params3[2] * 2. * params3[5]:
                    avera = params3[2]
            elif ng == 2:
                if params2[2] * 2. * params2[5]:
                    avera = params2[2]
        else:
            dAz = GalacticPolRa - avera
            if rs.telaz > 90. and rs.telaz < 270.:
                dEl = founddec - rs.dec
            else:
                dEl = rs.dec - founddec
        # 
        print("Ave RA: %7.3fd (%7.3fh)" % 
                  (avera, avera/15.))
        print("dAz: %7.3fd, dEl: %7.3fd" % 
              (dAz, dEl))
        
        writeDazDel( dazdelfilename, utca, utcb, lastId,
                     rs.telaz, rs.telel, dAz, dEl)

        aveutc, duration = radioastronomy.aveutcs(utc1, utc2)
        print("Average of crossing Times: %s, Time Interval: %8.2fs" % \
              (aveutc, duration))
        plt.plot(secs, tSums, color='blue',lw=3,
                 label='Intensities')
        if params2[0] != 0.:
            plt.plot(secs, bimodal(secs,*params2),color='red',lw=3,
                 label='2 Gaussians')
        if params3[0] != 0.:
            plt.plot(secs, trimodal(secs,*params3),color='green',lw=3,
                 label='3 Gaussians')
        if params4[0] != 0.:
            plt.plot(secs, quadmodal(secs,*params4),color='gold',lw=3,
                 label='4 Gaussians')
        if params5[0] != 0.:
            plt.plot(secs, fivemodal(secs,*params5),color='cyan',lw=3,
                 label='5 Gaussians')

        plt.xlabel( "Time (Seconds since Midnight %s)   dAz:%5.1fd dEl:%5.1fd" % (firstdate, dAz, dEl))
        plt.ylabel( "Integrated Intensity (K km/sec)")
        plt.title("%s Galactic Integrated Intensities - Tel:%2d Az:%7.1fd El:%7.1fd" % (
            firstdate, lastId, lastaz, lastel))
        plt.legend()

        # use the covariance matrix to get an estimate of fit uncertainty
        print("     Time    +/-   Intensity    +/-     Width   +/- ")
        i = 0
        if params2[0] != 0. and ng > 1:
            print("2 Gaussian Fit:")
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % ( 
                params2[0+i*3], sigma2[0+i*3], 
                params2[2+i*3], sigma2[2+i*3], 
                params2[1+i*3], sigma2[1+i*3]))
            i = 1
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % (
                params2[0+i*3], sigma2[0+i*3], 
                params2[2+i*3], sigma2[2+i*3], 
                params2[1+i*3], sigma2[1+i*3]))
            print("Delta: %8.1f  %7.1f (s)" % (
                params2[0]-params2[3],
                np.sqrt(sigma2[0]**2 + sigma2[3]**2)))
        i = 0
        if params3[0] != 0. and ng > 2:
            print("3 Gaussian Fit:")
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % ( 
                params3[0+i*3], sigma3[0+i*3], 
                params3[2+i*3], sigma3[2+i*3], 
                params3[1+i*3], sigma3[1+i*3]))
            i = 1
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % ( 
                params3[0+i*3], sigma3[0+i*3], 
                params3[2+i*3], sigma3[2+i*3], 
                params3[1+i*3], sigma3[1+i*3]))
            print("Delta: %8.1f  %7.1f (s)" % (
                params3[0]-params3[3], np.sqrt(sigma4[0]**2 + sigma4[3]**2)))
        i = 0
        if params4[0] != 0. and ng > 3:
            print("4 Gaussian Fit:")
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % ( 
                params4[0+i*3], sigma4[0+i*3], 
                params4[2+i*3], sigma4[2+i*3], 
                params4[1+i*3], sigma4[1+i*3]))
            i = 1
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % ( 
                params4[0+i*3], sigma4[0+i*3], 
                params4[2+i*3], sigma4[2+i*3], 
                params4[1+i*3], sigma4[1+i*3]))
            print("Delta: %8.1f  %7.1f (s)" % (
                params4[0]-params4[3], np.sqrt(sigma4[0]**2 + sigma4[3]**2)))
        i = 0
        if params5[0] != 0. and ng > 4:
            print("5 Gaussian Fit:")
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % ( 
                params5[0+i*3], sigma5[0+i*3], 
                params5[2+i*3], sigma5[2+i*3], 
                params5[1+i*3], sigma5[1+i*3]))
            i = 1
            print(" %8.1f  %7.1f   %7.1f  %7.1f %7.1f  %7.1f" % ( 
                params5[0+i*3], sigma5[0+i*3], 
                params5[2+i*3], sigma5[2+i*3], 
                params5[1+i*3], sigma5[1+i*3]))
            print("Delta: %8.1f  %7.1f (s)" % (
                params5[0]-params5[3], np.sqrt(sigma5[0]**2 + sigma5[3]**2)))
    
        plt.show()
        print( "Using %d values from file %s" % ( nData, filename))
           
    return
           

if __name__ == "__main__":
    main()
