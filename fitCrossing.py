"""
Read in an observation summary and fit the times of galaxy crossings.
From these measurements estimate the Azimuth and Elevation of the
telescope pointing.  Then compute the azimuth and elevation offsets
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 22May06 GIL iteratively fit peaks
# 22May04 GIL get azimuth, elevation and telescope position
# 22May02 GIL initial version

import sys
import os
import numpy as np
#from numba import jit
from matplotlib import pyplot as plt
from pylab import *
from scipy.optimize import curve_fit
import datetime
import gainfactor as gf

import radioastronomy

# ine file paths
offsetfilename = "/home/glen/GalacticPlaneCrossingOffsets.txt"
dazdelfilename = "/home/glen//2022-DazDel.txt"
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
    """
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

def writeDazDel( dazdelfilename, utc, cpuIndex, az, el, dAz, dEl):
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
    outline =  "%s %2d %8.2f %8.2f %6.2f %6.2f \r\n" % \
        (utc, cpuIndex, az, el, dAz, dEl)

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

    print( "dRa = %7.1f coresponds to dec %7.1fd (%d) in array %7.1f" %
           (dRa, decs[idec], idec, dRas[idec]))
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

    firstdate, utcs, tSums, dTs = gf.readAllValues( filename)

    # create an array of seconds
    nt = len(tSums)
    # now fit two gausians

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
                width = utcs[iMax] - utcs[j]
                break
            elif tTemps[k] < tHalf:
                width = utcs[k] - utcs[iMax]
                break
        # end searching for width
                
        print("t%d = %8.1f, y = %8.3f, width = %7.2f " % \
          (ng+1, utcs[iMax], tMax, width))
                
        # estimate the fit the gaussian
        expected=(utcs[iMax], width, tMax)

        sigma1 = [ 0., 0., 0.]
        try:
            params1,cov1=curve_fit(gauss,utcs,tTemps,expected,sigma=dTs)
            sigma1=sqrt(diag(cov1))
        except:
            print("Error fitting gaussian %d" % (ng+1))
            params1 = expected

        # fit was successful, subtract fit and try again
        for i in range(nt):
            tTemps[i] = tTemps[i] - gauss( utcs[i], *params1)

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
    tmin = utcs[0] - 10.
    tmax = utcs[nt-1] + 10.

    # now try 1, 2, 3 and 4 gaussians
    # keep the largest number that fits

    expected1 = ( utcPeaks[0], widths[0], tPeaks[0])
    sigma1 = [0., 0., 0.]
    try:
        params1,cov1=curve_fit(gauss,utcs,tSums,expected1,sigma=dTs)
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
        params2,cov2=curve_fit(bimodal,utcs,tSums,expected2,sigma=dTs,
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
        
        params3,cov3=curve_fit(trimodal,utcs,tSums,expected3,sigma=dTs,
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
        params4,cov4=curve_fit(quadmodal,utcs,tSums,expected4,sigma=dTs,
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
        params5,cov5=curve_fit(fivemodal,utcs,tSums,expected5,sigma=dTs,
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
    
    return firstdate, utcs, tSums, dTs, ng, params1, sigma1, params2, sigma2, params3, sigma3, params4, sigma4, params5, sigma5

def main():
    """
    Main executable for gridding astronomical data
    """

    nargs = len(sys.argv)

    if nargs < 2:
        print('Fit Crossings: find the time of crossings of peak intensities')
        print('fitCrossing <Spectrum Summary-files>')
        sys.exit()
    print( "reading %d files" % (nargs))

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
        firstdate, utcs, tSums, dTs, ng, params1, sigma1, params2, sigma2, params3, sigma3, params4, sigma4, params5, sigma5 = fitCrossing( filename)

        print("Fit successful for %d gausians" % (ng))
        lastId, lastaz, lastel = gf.lastazel()
        rs.telaz = lastaz
        rs.telel = lastel
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
        utc1 = utcmidnight
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
        
        utc0 = utcmidnight + datetime.timedelta( seconds = utcs[0])
        writeDazDel( dazdelfilename, utc0, lastId,
                     rs.telaz, rs.telel, dAz, dEl)

        aveutc, duration = radioastronomy.aveutcs(utc1, utc2)
        print("Average of crossing Times: %s, Time Interval: %8.2fs" % \
              (aveutc, duration))
        plt.plot(utcs, tSums, color='blue',lw=3,
                 label='Intensities')
        if params2[0] != 0.:
            plt.plot(utcs, bimodal(utcs,*params2),color='red',lw=3,
                 label='2 Gaussians')
        if params3[0] != 0.:
            plt.plot(utcs, trimodal(utcs,*params3),color='green',lw=3,
                 label='3 Gaussians')
        if params4[0] != 0.:
            plt.plot(utcs, quadmodal(utcs,*params4),color='gold',lw=3,
                 label='4 Gaussians')
        if params5[0] != 0.:
            plt.plot(utcs, fivemodal(utcs,*params5),color='orange',lw=3,
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
            print(" %8.1f  %7.1f" % ( 
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
            print(" %8.1f  %7.1f" % (
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
            print(" %8.1f  %7.1f" % (
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
            print(" %8.1f  %7.1f" % (
                params5[0]-params5[3], np.sqrt(sigma5[0]**2 + sigma5[3]**2)))
    
        plt.show()
        print( "Read %d values from file %s" % (len(utcs), filename))
           
    return
           

if __name__ == "__main__":
    main()
