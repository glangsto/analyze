#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#18FEB06 GIL revise for new name format, read date from file
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import matplotlib.pyplot as plt
import sys
import statistics
import radioastronomy
import numpy as np
from bandpass import *
from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter
dy = -1.

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

#for symbol, value in locals().items():
#    print symbol, value

nplot = 0
for iii in range(1, min(nargs,20)):

    filename = sys.argv[iii]

    rs = radioastronomy.Spectrum()
#    print filename
    rs.read_spec_ast( filename)
    rs.azel2radec()    # compute ra,dec from az,el

#    print("GAL Lon,Lat: %8.3f, %8.3f"  % (rs.gallon, rs.gallat))


    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    parts = aname.split('T')
    date  = parts[0]
    time  = parts[1]
    time  = time.replace('_',':')
    
    gallon = rs.gallon
    gallat = rs.gallat
    label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( time,rs.telaz,rs.telel,gallon,gallat)
    xv = rs.xdata  * 1.E-6 # convert to MHz
    yv = rs.ydataA * scalefactor
    # flag out the edges and middle
    nData = len(xv)
    nData2 = nData/2
# assume band edges match.  
    y0 = (yv[2]+yv[3]+yv[nData-3]+yv[nData-2])*.25
# replace edges with average
    yv[0] = y0
    yv[1] = y0
    yv[nData-1] = y0
# Now fix the middle down spike
    yv[nData2] = (yv[nData2-2]+yv[nData2+2])*.5
    yv[nData2-1] = ((3.*yv[nData2-2])+yv[nData2+2])*.25
    yv[nData2+1] = ((yv[nData2-2])+(3.*yv[nData2+2]))*.25

# now remove nominal zero
    yv = yv - y0

    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)
    ymin = min(yv)
    ymax = max(yv)
    ymed = statistics.median(yv)
    count = rs.count

    print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
    if nplot <= 0:
        fig,ax1 = plt.subplots(figsize=(10,6))
        plt.hold(True)
        fig.canvas.set_window_title(date)
        nplot = nplot + 1
    note = rs.noteA
#    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
    plt.ylim(0.9*ymin,1.1*yallmax)

    nData = len(xv)
    begin = int(nData2)
    end = nData-begin
#    p0 = estimatefit( xfit, yfit)
    spl = UnivariateSpline(xv, yv, s=100.)
    yfit = savgol_filter( yv, 101, 2)
    spl.set_smoothing_factor(.1)

#    yfit = spl(xv)
#
#    popt, pcov = curve_fit(bandpass, xfit, yfit, p0=p0)

#    yfit = bandpassparms(xv, popt)
    plt.plot( xv, yfit, 'b', linestyle='-', label='Smooth')

    plt.plot(xv, yv, colors[nplot], linestyle=linestyles[nplot],label=label)
#    plt.plot(xv, y_new, colors[nplot+1], linestyle=linestyles[nplot],label=label)
    nplot = nplot + 1
    
plt.title(note)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
