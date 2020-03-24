#Python Script to plot raw NSF events 
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#20Mar24 GIL increment color in plot order
#20Feb08 GIL date to title
#19NOV08 GIL add title option, switch to micro-seconds
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

linestyles = ['-','-.','--','-.','--','-','--','--','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9
mytitle = ""      # define plot title
doDebug = False
doAzEl = False
doMag = False

linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

if nargs < 2:
    print("E: Plot Events")
    print("Usage: ")
    print(" E [-A] [-T 'My Plot Title'] fileNames ")
    print("where:")
    print("  -T <plot title String> Enables user labeling of the plot")
    print("")
    print("Where many parameters are optional:")
    print("-A Show Az,El instead of Ra,Dec")
    print("-M Compute Magnitude")
    print("")
    print("Glen Langston - 2020 March 24")
    exit()

namearg = 1
iarg = 1          # start searching for input flags

# for all arguments, read list and exit when no flag argument found
while iarg < nargs:
    if sys.argv[iarg].upper() == '-T':   # now look for flags with arguments
        iarg = iarg+1
        mytitle = sys.argv[iarg]
        print('Plot title: ', mytitle)
    elif sys.argv[iarg].upper() == '-D':   # if debugging
        doDebug = True
        print("-D Additional Debug printing.")
    elif sys.argv[iarg].upper() == '-A':   # if labeling AzEl
        doAzEl = True
    elif sys.argv[iarg].upper() == '-M':   # if labeling AzEl
        doMag = True
    else:
        break
    iarg = iarg+1
namearg = iarg

firstel = -100.

note = "Events"
nplot = 0
for iii in range(namearg, min(nargs,25)):

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

    timestr = "%s" % (rs.utc)
    eMjd = rs.emjd
    iseconds = int(eMjd)
    seconds = 86400. * (eMjd - iseconds)

    if doDebug:
        print("%s == %s == %12.6f" %( time, timestr, seconds))
    parts = timestr.split('T')
    if len(parts) > 1:
        atime = parts[1]
    else:
        atime = timestr

    if doDebug:
        print("%s == %s" %( time, atime))

    gallon = rs.gallon
    gallat = rs.gallat
    ra = rs.ra
    dec = rs.dec
    if doAzEl:
        label = '%s A,E: %5.1f,%5.1f' % ( time,rs.telaz,rs.telel)
    else:
        label = '%s R,D: %6.2f,%6.2f, Lon,Lat=%5.1f,%5.1f' % ( time,rs.ra,rs.dec,gallon,gallat)
    xs = rs.xdata * 1.E6
    ya = rs.ydataA
    yb = rs.ydataB
    if rs.nTime < 1:
        print("Not an Event: ",filename)
        continue

    j = 0
    # convert time offsets to micro-seconds
    dt = 1.E6 * 0.5/rs.bandwidthHz
    if doMag:   # if plotting magnitude
        xv = np.zeros(rs.nSamples)
        yv = np.zeros(rs.nSamples)
        t = -dt * rs.refSample
        for i in range(rs.nSamples):
            y  = (ya[i]*ya[i]) + (yb[i]*yb[i])
            yv[j] = np.sqrt(y)
            xv[j] = t
            j = j + 1
            t = t + dt
    else:
        xv = np.zeros(rs.nSamples*2)
        yv = np.zeros(rs.nSamples*2)
        t = -2. * dt * rs.refSample
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

    n2 = int( rs.nSamples/2.5)
    ymin = min(yv)
    ymax = max(yv)
    yrms = np.std( yv[0:n2])

    count = rs.count

    if mytitle == "":
        mytitle = rs.site

    ypeak = max( -ymin, ymax)
    if yrms > 0:
        snr = ypeak/yrms
    else:
        snr = 0.

    print((' Ra: %6.2f Dec: %6.2f Max: %8.3f +/- %7.3f SNR: %6.1f ; %s %s' % (ra, dec, ypeak, yrms, snr, count, label)))
    if nplot <= 0:
        title = mytitle + " Az,El: %6.1f,%6.1f - %s" % (rs.telaz, rs.telel, date)
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(title)
    nplot = nplot + 1
    note = rs.noteA
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
    plt.ylim(yallmin,1.25*yallmax)


    plt.plot(xv, yv, colors[nplot-1], linestyle=linestyles[nplot-1],label=label)
plt.title(title)
plt.xlabel('Time Offset From Event (micro-seconds)')
plt.ylabel('Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
