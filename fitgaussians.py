"""
Read Astronomy files and fit a number of Gaussians to spectra
HISTORY
23Sep01 GIL return's x,y and width in FWHM units
23Aug31 GIL initial version

Glen Langston, National Science Foundation

"""
from scipy.optimize import curve_fit
import numpy as np

fwhmFactor = 1./(2*np.sqrt(np.log(2)))
# define the function of N gaussians

def func(x, *params):
    """
    This function computes the guassian values for an arbitrary
    number of gaussians.   The gaussians are assumed to have 3 parameters
    ctr = Center coordinate in units of x, the x axis fo the functin
    amp = Amplitude of the gaussian
    wid = Gaussian width, fwhmFactor converts to FWHM width values
    """
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = fwhmFactor*params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

def fitgaussians( Ngauss, x, y, minSigma = 2, \
                  doFinal = True, doPlot = False, verbose = False):
    """
    fitgaussians() fits up to Ngaussians and returns fit parameters
    The fitting process stops if fitting fails or the peak is smaller
    than the noise.   The fitting is done interatively, finding and subtracting
    fits until the remaining peak is less than RMS*minSigma
    Where 
    Ngauss      - maximm number of gaussians to fit
    x,y         - arrays of x and y values to fit
    minSigma    - stop trying to find gaussians if peak less than sigma
    doFinal     - fit final merged collection of gaussians.
                  with a large number of gaussians this yields strange results
    doPlot      - optionally plot steps of fitting
    verbose     - optionally print lots of debugging info.
    """

    if doPlot:
        try:
            import matplotlib.pyplot as plt
        except:
            doPlot = False

    # estimate the fit width as a fraction of the x range
    width = np.max(x) - np.min(x)
    width = width/(20.* Ngauss)

    # copy y data for finding guesses
    yguess = y
    nfit = 0
    for ng in range( Ngauss):
        # find min and max to use one as guess
        iymax = np.argmax(yguess)
        iymin = np.argmin(yguess)
        ymax = yguess[iymax]
        xmax = x[iymax]
        ymin = yguess[iymin]
        xmin = x[iymin]
        rms = np.std(yguess)

        # test whether fitting should stop
        peak = ymax
        if peak < - ymin:
            peak = -ymin
        # if time to exit, no peaks left
        if peak < minSigma * rms:
            if verbose:
                print( "Fitting stopped, remaining peak is small")
                print( "Peak %.2f <  %.1f * %.2f" % (peak, minSigma, rms))
            break
                       
        if verbose:
            print( "Min %.2f at %.2f;  Max %.2f at %.2f; RMS = %.2f" % \
                   (ymin, xmin, ymax, xmax, rms))

        # fit the negative peak, if bigger
        if ymax > -ymin : 
            guess = [xmax, ymax, width]
        else:
            guess = [xmin, ymin, width]

        # now fit this gaussian to the residual 
        try:
            popt, pcov = curve_fit(func, x, yguess, p0=guess)
            nfit = ng + 1
            # create an array uncertaintys
            rmsfit = np.zeros(3)
            rmsfit[0] = np.sqrt(pcov[0,0])
            rmsfit[1] = np.sqrt(pcov[1,1])
            rmsfit[2] = np.sqrt(pcov[2,2])
        except:
            print("Fitting %d Gaussian failed, exiting!" % (ng+1))
            nfit = ng
            break

        # now transfer only the fit to the combined input
        if ng == 0:
            allfit = np.zeros(3*Ngauss)
            allrms = np.zeros(3*Ngauss)
        # now transfer results
        ng3 = 3 * ng
        # now transfere results and rms
        for iii in range(3):
            allfit[ng3+iii] = popt[iii]
            allrms[ng3+iii] = rmsfit[iii]

        if verbose:
            print("Fit of %d gaussians: " % (nfit))
            print(popt)
            print("+/-")
            print(rmsfit)
            
        fit = func(x, *popt)
        if doPlot:
            plt.plot(x, y)
            plt.plot(x, yguess)
            plt.plot(x, fit)
            plt.show()
        # next cycle fit one more gaussians to residual
        # guess finding is done by subtracting preceding fit
        yguess = yguess - fit
        # end of finding guesses to fit
        
    if nfit > 0 and doFinal:
        print("Performing simultaneous fit of %d gaussians" % (nfit))
        
        # now trim to number of items fit
        allfit = allfit[0:nfit*3]
        allrms = allrms[0:nfit*3]
        # now guessing is complete; fit original data with the guess
        try:
            # in final fit, work with original y values
            popt, pcov = curve_fit(func, x, y, p0=allfit, \
                                   maxfev=3000 + nfit*1000)
            if verbose:
                print("Simultaneous Fit Successful")
            fitOK = True
            for iii in range(nfit*3):
                allfit[iii] = popt[iii]
                allrms[iii] = np.sqrt(pcov[iii,iii])
        except:
            if verbose:
                print("Combined Curve_Fit failed, returning sum of fits")
            fitOK = False

    # now trim to number of reasonable gaussians fit
    allfit = allfit[0:nfit*3]
    allrms = allrms[0:nfit*3]
    if doPlot:
        # now comput final fit
        fit = func(x, *allfit)
        plt.plot(x, y)
        plt.plot(x, fit , 'r-')
        plt.show()
    # return the number of gaussians fit and parameters
    return nfit, allfit, allrms
    # end of fitgaussians()
    
def main():
    """
    Test executable 
    Takes as input the number of gaussians and a list of files to fit
    """
    import sys
    import numpy as np
    import radioastronomy
    import gainfactor
    doPlot = True
    
    nargs = len(sys.argv)
    if nargs < 2:
        print('fitgaussians: N-Gaussians N-Sigma, Filenames')
        print('where')
        print(' N-Gaussians - number of gaussians to fit to the spectra, ie 4')
        print(' N-Sigma     - limit on peak to fit, relative to noise, ie 3')
        print(' Filenames   - Radio Astronomy Files to fit gaussians')
        print(' ')
        print('Fitting is better if a baseline is subtracted from the data')
        print('')
        print('Glen Langston,  National Science Foundation, 2023 August 31')
        exit()

    # define min and max velocities for fitting (later these will be arguemnts)
    maxvel = 200000.
    minvel = -maxvel
    
    Ngauss = int(sys.argv[1])
    if Ngauss < 1:
        print('Number of gaussians must be greater than 0: %s' % \
              (sys.argv[0]))
        exit()
    minSigma = float(sys.argv[2])
    if minSigma < .1:
        print("N-sigma too small for reliable fits: %f" % (minSigma))
    names = sys.argv[3:] 

    # prepare to read the spectra
    rs = radioastronomy.Spectrum()

    for filename in names:
        rs.read_spec_ast(filename)
        rs.azel2radec()
        # kelvins files are usually intensity versus velocity
        vel = rs.xdata
        imin, imax = gainfactor.velocity_to_indicies( vel, minvel, maxvel)
        y = rs.ydataA[imin:imax]
        x = vel[imin:imax]/1000.   # convert to km/sec
        # guess width as input velocity range
        width = vel[imin] - vel[imax]

        nfit, allfit, allrms = fitgaussians( Ngauss, x, y,
                                             minSigma=minSigma, \
                                             doFinal=False, \
                                             doPlot=doPlot, verbose=False)
        
        for iii in range(nfit):
            iii3 = 3*iii
            print("Peak %2d: x,y = %7.1f+/-%7.1f, %7.1f+/-%7.1f, FWHM = %7.1f+/-%7.1f" % \
                  (iii+1, allfit[iii3], allrms[iii3], \
                   allfit[iii3+1], allrms[iii3+1],  \
                   allfit[iii3+2], allrms[iii3+2]))
        # end of fitting one file
    # end of main

if __name__ == "__main__":
    main()
