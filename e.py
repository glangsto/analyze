#Python Script to plot raw NSF events 
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#23Feb24 GIL keep track of min and max
#23Feb17 GIL fix np deletion of np.int and np.float
#22Apr07 GIL clean up messages
#22Mar23 GIL pull out leading pi from file name
#21Aug21 GIL deal with different length event files
#21Aug05 GIL matplotlib updates
#21Mar21 GIL initialize t variable for computing M
#20Dec10 GIL optionally shift plots in y 
#20Aug28 GIL optionally on plot to files
#20Mar24 GIL increment coyalllor in plot order
#20Feb08 GIL date to title
#19NOV08 GIL add title option, switch to micro-seconds
#19OCT10 GIL also show ra and dec
#19MAR11 GIL fix plot axies
#19JAN16 GIL initial version
#
import matplotlib as mpl
#import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
import numpy as np

# set default offset
dy = -1.
# set default telescope index
intel = 0
nargs = len( sys.argv)

linestyles = ['-','-.','-', '-', '--','-.','--','-','--','--','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['b','r','g','c', 'm', 'b','r','g','b','c', 'm', 'r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g']

xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9
mytitle = ""      # define plot title
doDebug = False
doAzEl = False
doMag = False
fileTag = ""
yoffset = 0
y0 = 0

firstdate = ""
xa = -1
xb = -1
lastn = -1

if nargs < 2:
    print("E: Plot Events")
    print("Usage: ")
    print(" E [-A] [-P] [-T 'My Plot Title'] fileNames ")
    print("where:")
    print("  -T <plot title String> Add Title to the plot")
    print("  -P write PNG and PDF files instead of showing plot")
    print("  -A Show Az,El instead of Ra,Dec")
    print("  -B <sample> Set first sample to plot (default is 1/4 of samples)")
    print("  -I <telescope index> Set the telescope identifier number")
    print("  -E <sample> Set last sample to plot (default is end of samples)")
    print("  -M Compute Magnitude")
    print("  -Y <offset> optionally add an offset to each plot")
    print("  -Z <file tag> optionally add tag to PDF and PNG file names")
    print("")
    print("Glen Langston - NSF    2020 August 28")
    exit()

namearg = 1
iarg = 1          # start searching for input flags
doPlotFile = False

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
    elif sys.argv[iarg].upper() == '-B':   # if setting beginning sample
        iarg = iarg + 1
        xa = int( sys.argv[iarg])
    elif sys.argv[iarg].upper() == '-E':   # if setting ending sample
        iarg = iarg + 1
        xb = int( sys.argv[iarg])
    elif sys.argv[iarg].upper() == '-I':   # if setting ending sample
        iarg = iarg + 1
        intel = int( sys.argv[iarg])
    elif sys.argv[iarg].upper() == '-M':   # if labeling AzEl
        doMag = True
    elif sys.argv[iarg].upper() == '-P':
        doPlotFile = True
    elif sys.argv[iarg].upper() == '-Y':   # if offsetting in Y
        iarg = iarg + 1
        yoffset = float( sys.argv[iarg])
    elif sys.argv[iarg].upper() == '-Z':     # label written files
        iarg = iarg+1
        fileTag = str(sys.argv[iarg])
        print( 'File tag: %s' % (fileTag))
    else:
        break
    iarg = iarg+1


# end of reading arguments
# to create plots in cronjobs, must use a different backend
if doPlotFile:
    mpl.use('Agg')
import matplotlib.pyplot as plt
import gainfactor as gf

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

# if the pi??-events name is present get index
    itel = intel
    piparts = filename.split('-events')
    # if '-events' in the file 
    if len(piparts) > 0:
        telparts = piparts[0]
#        print(telparts)
        # if directory starts with pi
        if telparts[0:2] == "pi":
            # keep the integer value
            telparts = telparts[2:]
#            print(telparts)
            itel = int(telparts)
            
    if firstdate == "":
        firstdate = date
    lastdate = date
    
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
    if itel > 0:
        label = ("%d: " % (itel)) + label
    if (xa < 0) or (lastn != rs.nSamples):
          xa = int(rs.nSamples/4)
    if (xb < 0) or (lastn != rs.nSamples):
          xb = rs.nSamples-1
    lastn = rs.nSamples
    
    xs = rs.xdata[xa:xb] * 1.E6

    ya = rs.ydataA[xa:xb]
    yb = rs.ydataB[xa:xb]
    nData = xb-xa
    if rs.nTime < 1:
        print("Not an Event: ",filename)
        continue

    j = 0
    # convert time offsets to micro-seconds
    dt = 1.E6 * 0.5/rs.bandwidthHz
    # compensate for starting earlier in plot
    if doMag:   # if plotting magnitude
        xv = np.zeros(nData)
        yv = np.zeros(nData)
        t = (2 * -dt * (rs.refSample-xa))
        for i in range(nData):
            y  = (ya[i]*ya[i]) + (yb[i]*yb[i])
            yv[j] = np.sqrt(y)
            xv[j] = t
            j = j + 1
            t = t + (2.*dt)
    else:
        xv = np.zeros(nData*2)
        yv = np.zeros(nData*2)
        t = (-2. * dt * (rs.refSample-xa))
        for i in range(nData):
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

    n2 = int( nData/2.5)
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
#        fig.set_window_title(title)
    nplot = nplot + 1
    note = rs.noteA
    yallmin = min(ymin+y0,yallmin)
    yallmax = max(ymax+y0,yallmax)

    plt.plot(xv, yv+y0, colors[nplot-1], linestyle=linestyles[nplot-1],label=label)
    y0 = y0 + yoffset
plt.title(title)
plt.xlabel('Time Offset From Event (micro-seconds)')
plt.ylabel('Intensity (Counts)')
plt.legend(loc='upper right')
plt.xlim(xallmin,xallmax)
if not doMag:
#    if y0 > 0.:
#        plt.ylim(yallmin,yallmax+y0*.666)
#    elif y0 < 0.:
#        plt.ylim(yallmin+y0*.666,yallmax)
#    else:
    plt.ylim(yallmin,(1.1*yallmax))
if doPlotFile:
    if fileTag == "":
        fileTag = "E-" + firstdate
    outpng = fileTag + ".png"
    plt.savefig(outpng,bbox_inches='tight')
    outpdf = fileTag + ".pdf"
    plt.savefig(outpdf,bbox_inches='tight')
    print( "Wrote files %s and %s" % (outpng, outpdf))
else:
    # else show the plots
    plt.show()
