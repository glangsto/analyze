"""
Read in an observation summary and fit the times of galaxy crossings.
From these measurements estimate the Azimuth and Elevation of the
telescope pointing.  Then compute the azimuth and elevation offsets
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 22May04 GIL get azimuth, elevation and telescope position
# 22May02 GIL initial version

import sys
import copy
import numpy as np
#from numba import jit
from matplotlib import pyplot as plt
from pylab import *
from scipy.optimize import curve_fit
import datetime
import gainfactor as gf

#from matplotlib import colors
import radioastronomy

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

    params2,cov2=curve_fit(bimodal,utcs,tSums,expected2)
    sigma2=sqrt(diag(cov2))

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

    params3,cov3=curve_fit(trimodal,utcs,tSums,expected3)

    params4,cov4=curve_fit(quadmodal,utcs,tSums,expected4)
    
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

    # create a radio astronomy spectrum structure to use utilities
    rs = radioastronomy.Spectrum()

# first read through all data and find hot load
    names = sys.argv[1:]
    names = sorted(names)

    print(names)
    firsttime = ""
    lasttime = ""
    firstdate = ""
    lastdate = ""
    date = "foo"
    count = 0
    utcs = []

 #   print( "Reading %d files" % (len(names)))
           
    for filename in names:
        print ("File: %s" % (filename))

        # return the data and three fits
        firstdate, utcs, tSums, dTs, \
            params2, cov2, params3, cov3, params4, cov4 = \
            fitCrossing( filename)
    
        plt.plot(utcs, tSums, color='blue',lw=3,
                 label='Intensities')
        plt.plot(utcs, bimodal(utcs,*params2),color='red',lw=3,
                 label='2 Gaussians')
        plt.plot(utcs, trimodal(utcs,*params3),color='green',lw=3,
                 label='3 Gaussians')

        plt.xlabel( "Time of Day (Seconds since Midnight %s)" % (firstdate))
        plt.ylabel( "Integrated Intensity (K km/sec)")
        lastId, lastaz, lastel = gf.lastazel()
        plt.title("%s Galactic Integrated Intensities - Tel:%2d Az:%7.1fd El:%7.1fd" % (
            firstdate, lastId, lastaz, lastel))
        plt.legend()
        plt.show()

        # use the covariance matrix to get an estimate of fit uncertainty
        sigma2=sqrt(diag(cov2))
        sigma3=sqrt(diag(cov3))
        #sigma4=sqrt(diag(cov4))
        
        print("Estimated time interval between Galactic Plane Crossings")
        print("2 Gaussian Fit:")
        print("     Time          Intensity            Width       ")
        i = 0
        print(" %8.1f  %5.1f   %7.1f  %5.1f  %7.2f   %5.2f" % ( 
            params2[0+i*3], sigma2[0+i*3], 
            params2[1+i*3], sigma2[1+i*3], 
            params2[2+i*3], sigma2[2+i*3]))
        i = 1
        print(" %8.1f  %5.1f   %7.1f  %5.1f  %7.2f  %5.2f" % (
            params2[0+i*3], sigma2[0+i*3], 
            params2[1+i*3], sigma2[1+i*3], 
            params2[2+i*3], sigma2[2+i*3]))
        print(" %8.1f  %5.1f" % ( 
            params2[0]-params2[3], np.sqrt(sigma2[0]**2 + sigma2[3]**2)))
        print("3 Gaussian Fit:")
        i = 0
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
