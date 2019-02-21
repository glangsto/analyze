#python program to take two arrays and find the RFI
import numpy as np
import statistics

from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def gaus(x,a,x0,sigma):
    """
    gaus() computes the gausian value at x, assuming the  model input values
    inputs:
    x     coordinate for which the gaussian value is computed
    a     gausian model height
    x0    coordinates of the center of the gaussian
    sigma width of the gaussian model
    """    
    return a*exp(-(x-x0)**2/(2*sigma**2))

def findRfi(xv, yv, nMedian=17):
    """
    findRfi(xv, yv, nMedian) finds RFI and separates the narrow lines from the median of the spectra
    findRfi takes two input arrays and finds the RFI in the array, assuming the RFI bandwidth is only a few channels
    findRfi returns two arrays, which should equal the input array plus the remaining signal
    """

    nData = len(yv)
    hv = np.zeros(nData)
    ys = np.zeros(nData)
    xs = np.zeros(nData)
    medianValues = np.zeros(nData)
    dyd = np.zeros(nData)
    myv = np.zeros(nData)
    nloop = 0
    
    for iii in range(nData):
        hv[iii] = yv[iii]

# loop through hot load removing RFI and writing RFI list
    while nMedian > 7:

#    print 'N median: %d ' % (nMedian)
        for iii in range(nMedian,nData-(nMedian-1)):
            ys = hv[(iii-nMedian):(iii+nMedian)]
            myv[iii] = statistics.median(ys)

        for iii in range(nMedian+1):
            myv[iii] = myv[nMedian+2]
            myv[nData-(iii+1)] = myv[nData-(nMedian+3)]

        for iii in range(nData):
            medianValues[iii] = medianValues[iii] + myv[iii]

        hrms = statistics.stdev(hv)
        # now subtract the median values
        hv = hv - myv

    # flatten out the ends

        dyrms = statistics.stdev(hv)

        print "RMS before and after median: ", hrms, dyrms, nMedian

        for iii in range(nMedian,int(nData-(nMedian+1)), int(nMedian/2)):
            # if a peak was found
            if (hv[iii] > 1.*dyrms) or (hv[iii] < -1.*dyrms):
#            print 'range: %d to %d (%d) ' % (iii-nMedian, iii+nMedian, nMedian)
                # select out a range of x and y values for a gaussian fit
                xs = xv[(iii-nMedian):(iii+nMedian)]
                ys = hv[(iii-nMedian):(iii+nMedian)]
                n = len(xs)
                x0 = xv[iii]
                y0 = hv[iii]
#                print 'Array sizes: ', n, y0, x0, (xv[iii]-xv[iii-2])
# now fit the gaussian                
                popt,pcov = curve_fit(gaus,xs,ys,p0=[y0, x0, (xv[iii]-xv[iii-2])])
                print iii, popt, '+/-', pcov[0]
# now subtract the fit
                for jjj in range( -nMedian,nMedian):
                    # compute the correction to the fit
                    ady = gaus( xv[jjj+iii], popt[0], popt[1], popt[2])
                    hv[jjj+iii] = hv[jjj+iii] - ady

# step through the median values
        nloop = nloop+1
        nMedian = nMedian - 2
    return hv, medianValues    

