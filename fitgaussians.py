"""
Read Astronomy files and fit Gaussians to spectra
"""
from scipy.optimize import curve_fit
import numpy as np

# define the function of N gaussians

def func(x, *params):
    """
    This function computes the guassian values for an arbitrary
    number of gaussians.   The gaussians are assumed to have 3 parameters
    ctr = Center coordinate in units of x, the x axis fo the functin
    amp = Amplitude of the gaussian
    wid = Gaussian width,  without the normalizing factor (2 pi)
    """
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

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
    try:
        import matplotlib.pyplot as plt
    except:
        doPlot = False
    verbose = False
    
    nargs = len(sys.argv)
    if nargs < 2:
        print('fitgaussians: N-Gaussians Filenames')
        print('where')
        print(' N-Gaussians - number of gaussians to fit to the spectra, ie 4')
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
    names = sys.argv[2:] 

    rs = radioastronomy.Spectrum()

    for filename in names:
        rs.read_spec_ast(filename)
        rs.azel2radec()
        # kelvins files are usually intensity versus velocity
        vel = rs.xdata
        imin, imax = gainfactor.velocity_to_indicies( vel, minvel, maxvel)
        y = rs.ydataA[imin:imax]
        x = vel[imin:imax]
        # guess width as input velocity range
        width = vel[imin] - vel[imax]
        if width < 0.:
            width = - width

        width = width/(20.* Ngauss)

        # copy y data for finding guesses
        yguess = y
        for ng in range( Ngauss):
            imax = np.argmax(yguess)
            ymax = y[imax]
            xmax = x[imax]
            if ng == 0:
                guess = [xmax, ymax, width]
            else:
                guess = guess + [xmax, ymax, width]
            # now fit this number of gaussians   

        
            try:
                popt, pcov = curve_fit(func, x, yguess, p0=guess)
                nfit = ng + 1
            except:
                print("Fitting %d gaussians failed!" % (ng+1))
                nfit = ng
                guess = guess[0:(3*nFit)]
                break
            print("Fit of %d gaussians: " % (nfit))
            fit = func(x, *popt)
            print(popt)
            if doPlot:
                plt.plot(x,yguess)
                plt.show()
            yguess = yguess - fit
            # end of finding guesses to fit
        
        if nfit > 0:
            print("Now performing final fit of %d gaussians" % (nfit))
        
            # now guessing is complete; fit original data with the guess
            popt, pcov = curve_fit(func, x, y, p0=guess)
            print("Final Fit:")
            print(popt)
            print("Final Covariance")
            print(pcov)
            if doPlot:
                plt.plot(x, y)
                plt.plot(x, fit , 'r-')
                plt.show()
        else:
            print("fitting gaussians failed!")
            
        # end of fitting one file
    # end of main

if __name__ == "__main__":
    main()
