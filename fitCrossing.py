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

# define file paths
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

def fitCrossing( filename):
    """
    fitCrossing takes the integrated intensities from a sumary file,
    then fits a series of gaussians
    """

    firstdate, utcs, tSums, dTs = gf.readAllValues( filename)

    # create an array of seconds
    nt = len(tSums)
    # now fit two gausians

    # get the indexs to  the first intensity peak
    imax1 = np.argmax( tSums)
    tMax1 = tSums[imax1]
    # now that one peak is found, must find 2nd peak far from first
    nt10 = int(nt/10)
    if imax1 > nt/2 :
        imax2 = np.argmax( tSums[0:(imax1-nt10)])
        tMax2 = tSums[imax2]
    else:
        imax2 = np.argmax( tSums[(imax1+nt10):nt])
        imax2 = imax2 + imax1 + nt10
        if imax2 >= nt:
            imax2 = nt - nt10
        tMax2 = tSums[imax2]
    
    print("t1 = %8.1f, y1 = %8.3f, width = %7.2f " % \
          (utcs[imax1], tMax1, 1800.))
    print("t2 = %8.1f, y2 = %8.3f, width = %7.2f " % \
          (utcs[imax2], tMax2, 1800.))
    print("Read %d samples from file: %s" % (nt, filename))
          
    #find the mid point intensity
    imax3 = int((imax1 + imax2)/2)
    tMax3 = tSums[imax3]

    if imax1 > imax2:
        imax4 = nt + imax1 - imax2
    else:
        imax4 = nt + imax2 - imax1
    if imax4 > nt:
        imax4 = imax4 - nt
    if imax4 < 0:
        imax4 = imax4 + nt
    tMax4 = tSums[imax4]

    # NOW pick the ihigher third point
    if tMax4 > tMax3:
        tMax3 = tMax4
        imax3 = imax4

    # first try 2 gaussians
    expected2=(utcs[imax1], tMax1, 1800., utcs[imax2], tMax2, 3600.)


    # now try 3 gaussians
    expected3=(utcs[imax1], tMax1, 1800., \
               utcs[imax2], tMax2, 1800., \
               utcs[imax3], tMax3, 3600.)

    # finally try 4 gaussians
    expected4=(utcs[imax1], tMax1, 1800., \
               utcs[imax2], tMax2, 1800., \
               utcs[imax3], tMax3, 3600., \
               utcs[imax3+1], tMax3/10, 7200.)
    
    fig, ax1 = plt.subplots(figsize=(10, 6))

    try:
        params3,cov3=curve_fit(trimodal,utcs,tSums,expected3,sigma=dTs)
    except:
        print("Error trying a 3 gaussian fit")
        params3 = [ 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        cov3 = params3
        
    if doTwo:
        try:
            params2,cov2=curve_fit(bimodal,utcs,tSums,expected2,sigma=dTs)
        except:
            print("Error trying a 2 gaussian fit")
            params2 = params3[0:6]
            cov2 = cov3[0:6]
    else:
        params2 = params3[0:6]
        cov2 = cov3[0:6]

    if doFour:
        try:
            params4,cov4=curve_fit(quadmodal,utcs,tSums,expected4, sigma=dTs)
        except:
            print("Error trying a 4 gaussian fit")
            params4 = [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
            cov4 = params4
    else:
        params4 = [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        cov4 = params4
    
    return firstdate, utcs, tSums, dTs, params2, cov2, params3, cov3, params4, cov4

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
        firstdate, utcs, tSums, dTs, \
            params2, cov2, params3, cov3, params4, cov4 = \
            fitCrossing( filename)

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
        if params4[0] != 0.:
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

#        if ra1 > 180.:
#            ra1 = ra1 - 360.
#        if ra2 > 180.:
#            ra2 = ra2 - 360.
        avera = (ra1 + ra2)/2.
        print("Ave RA: %7.3fd (%7.3fh)" % 
              (avera, avera/15.))
        dAz = GalacticPolRa - avera
        dRa = ra1 - ra2
        if dRa < 0:
            dRa = - dRa
        # read in offsets vs dec
        decs, ra1s, ra2s, dRas = readGalacticOffsets( offsetfilename)
        founddec = decFromCrossingDelta( dRa, decs, ra1s, ra2s, dRas)
        if rs.telaz > 90. and rs.telaz < 270.:
            dEl = founddec - rs.dec
        else:
            dEl = rs.dec - founddec
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
            plt.plot(utcs, quadmodal(utcs,*params4),color='green',lw=3,
                 label='4 Gaussians')

        plt.xlabel( "Time of Day (Seconds since Midnight %s)" % (firstdate))
        plt.ylabel( "Integrated Intensity (K km/sec)")
        plt.title("%s Galactic Integrated Intensities - Tel:%2d Az:%7.1fd El:%7.1fd" % (
            firstdate, lastId, lastaz, lastel))
        plt.legend()
        plt.show()

        # use the covariance matrix to get an estimate of fit uncertainty
        if params2[0] != 0.:
            sigma2=sqrt(diag(cov2))
        if params3[0] != 0.:
            sigma3=sqrt(diag(cov3))
        #sigma4=sqrt(diag(cov4))
        
        print("Estimated time interval between Galactic Plane Crossings")
        print("     Time          Intensity            Width       ")
        i = 0
        if params2[0] != 0.:
            print("2 Gaussian Fit:")
            print(" %8.1f  %5.1f   %7.1f  %5.1f  %7.2f  %5.2f" % ( 
                params2[0+i*3], sigma2[0+i*3], 
                params2[1+i*3], sigma2[1+i*3], 
                params2[2+i*3], sigma2[2+i*3]))
            i = 1
            print(" %8.1f  %5.1f   %7.1f  %5.1f  %7.2f  %5.2f" % (
                params2[0+i*3], sigma2[0+i*3], 
                params2[1+i*3], sigma2[1+i*3], 
                params2[2+i*3], sigma2[2+i*3]))
            print(" %8.1f  %5.1f" % ( 
                params2[0]-params2[3],
                np.sqrt(sigma2[0]**2 + sigma2[3]**2)))
        i = 0
        if params3[0] != 0.:
            print("3 Gaussian Fit:")
            print(" %8.1f  %5.1f   %7.1f  %5.1f  %7.2f  %5.2f" % ( 
                params3[0+i*3], sigma2[0+i*3], 
                params3[1+i*3], sigma2[1+i*3], 
                params3[2+i*3], sigma2[2+i*3]))
            i = 1
            print(" %8.1f  %5.1f   %7.1f  %5.1f  %7.2f  %5.2f" % ( 
                params3[0+i*3], sigma3[0+i*3], 
                params3[1+i*3], sigma3[1+i*3], 
                params3[2+i*3], sigma3[2+i*3]))
            print(" %8.1f  %5.1f" % (
                params3[0]-params3[3], np.sqrt(sigma3[0]**2 + sigma3[3]**2)))
    
        print( "Read %d values from file %s" % (len(utcs), filename))
           
    return
           

if __name__ == "__main__":
    main()
