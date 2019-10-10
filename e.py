#Python Script to plot raw NSF events 
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#19OCT10 GIL also show ra and dec
#19MAR11 GIL fix plot axies
#19JAN16 GIL initial version
#
import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
import numpy as np

dy = -1.

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

#for symbol, value in locals().items():
#    print symbol, value
firstel = -100.

note = "Events"
nplot = 0
for iii in range(1, min(nargs,25)):

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
#    time  = time.replace('_',':')
    
    gallon = rs.gallon
    gallat = rs.gallat
    ra = rs.ra
    dec = rs.dec
    label = '%s, RA,Dec: %7.3f,%7.3f, Lon,Lat=%5.1f,%5.1f' % ( time,rs.ra,rs.dec,gallon,gallat)
    xs = rs.xdata 
    ya = rs.ydataA
    yb = rs.ydataB
    if rs.nTime < 1:
        print "Not an Event: ",filename
        continue

    xv = np.zeros(rs.nSamples*2)
    yv = np.zeros(rs.nSamples*2)
    j = 0
    dt = 0.5/rs.bandwidthHz
    t = -2. * dt * rs.refSample
    if iii == -1:  # no op
        print "First Time %12.9f (s); Delta T = %12.9f (s)" % (t, dt)
    for i in range(rs.nSamples):
        yv[j] = ya[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        yv[j] = yb[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)

    ymin = min(yv)
    ymax = max(yv)
    ymed = (ymin+ymax)*.5
    count = rs.count

    print(' Ra: %8.5f Dec: %8.5f Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ra, dec, ymax, ymed, ymax/ymed, count, label))
    if nplot <= 0:
        title = date + " Az,El: %6.1f,%6.1f" % (rs.telaz, rs.telel)
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(title)
        nplot = nplot + 1
    note = rs.noteA
#    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
#    plt.ylim(0.9*ymin,1.5*yallmax)
#    plt.ylim(0.9*ymin,1.25*yallmax)
    plt.ylim(yallmin,1.25*yallmax)


    plt.plot(xv, yv, colors[iii-1], linestyle=linestyles[iii-1],label=label)
plt.title(note)
plt.xlabel('Time Offset From Event (s)')
plt.ylabel('Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
