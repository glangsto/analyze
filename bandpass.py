# script to simulate, then fit a SDR dongle band pass response
# HISTORY
# 16NOV15 GIL initial version
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def estimatefit( x, y):
     """
     Estimate fit returns 8 parameters for a model of band pass
     """
     n = len(x)
     n2 = n/2
     n3 = int(n/3)
     n4 = int(n/4)
     n8 = int(n/8)
     y0 = max(y[(n2-n8):(n2+n8)])
     y = y/y0
     imax = 1
     for iii in range(n8,n3):
          if y[iii] > y[imax]:
               imax = iii
     xa = x[imax]
     a = y[imax]
     # now find halfwidth
     for iii in range(n8,n3):
          if 2*y[iii] > y[imax]:
               imax = iii
               break
     dx = xa-x[imax]
     sa = 1./(dx*dx)
     # now start finding third peak
     imax = 2*n3
     for iii in range(2*n3,n-n8):
          if y[iii] > y[imax]:
               imax = iii
     xc = x[imax]
     c = y[imax]
     imax = n-n8
     for iii in range(n-n8,2*n3,-1):
          if 2.*y[iii] > y[imax]:
               imax = iii
               break
     dx = xc - x[imax]
     sc = 1./(dx*dx)
     # for middle peak just assume 1/8 of the width
     dx = x[n2]-x[n2-n8]
     sb = 1./(dx*dx)
     return xa, xc, a, y0, c, sa, sb, sc

def bandpassparms(x, parms):
     xa = parms[0]
     xc = parms[1]
     a = parms[2]
     y0 = parms[3]
     c = parms[4]
     sa = parms[5]
     sb = parms[6]
     sc = parms[7]
     return bandpass( x, xa, xc, a, y0, c, sa, sb, sc)

def bandpass( x, xa, xc, a, y0, c, sa, sb, sc):
     n = len(x)
     xb = x[int(n/2)]
     y = (a*np.exp(-sa*((x-xa)**2))) + (np.exp(-sb*((x-xb)**2))) + (c*np.exp(-sc*((x-xc)**2))) 
     return y*y0

def testbandpass():
     xa = 200. 
     xc = 1024-xa
     a = 1.3
     y0 = 1000.
     c = 1.3
     sa = 1./100.**2
     sb = 1./200.**2
     sc = 1./100.**2
     pinit = [xa, xc, a, y0, c, sa, sb, sc]

     xdata = np.linspace(0, 1024., 1024)

     y = bandpass(xdata, xa, xc, a, y0, c, sa, sb, sc)
     ydata = y + 100. * np.random.normal(size=len(xdata))

     plt.plot( xdata, ydata)
     p0 = estimatefit( xdata, ydata)
     popt, pcov = curve_fit(bandpass, xdata, ydata, p0=p0)

     print 'popt: ', popt
     print 'pcov: ', pcov

     xa = popt[0]
     xc = popt[1]
     a = popt[2]
     y0 = popt[3]
     c = popt[4]
     sa = popt[5]
     sb = popt[6]
     sc = popt[7]

     yfit = bandpass(xdata, xa, xc, a, y0, c, sa, sb, sc)
     plt.plot( xdata, yfit, 'r', linestyle='-')
     plt.show()

