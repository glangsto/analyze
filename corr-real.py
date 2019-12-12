#Python Script to correlate samples in raw NSF events 
#plot the self correlation values samples near an event
#HISTORY
#19DEC11 GIL test flipping the I/Q order with the -O argument
#19DEC10 GIL initial version
#
import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
import numpy as np
import pandas as pd

dy = -1.

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9
mytitle = ""      # define plot title
doAbs = False
doDebug = False
doZero = False
doOrder = False

linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

if nargs < 2:
    print "CORR: Cross correlate Event Samples over a range in sample offsets"
    print "Usage: "
    print " CORR [-T 'My Plot Title'] [-A -D -O -Z] <nSamples> fileNames "
    print "where:"
    print " <nSamples> is the maximum number of samples to compute the cross correlation"
    print "  A range of +/- nSamples correllations are performed and plotted"
    print "  -T <plot title String> Optionally Enables user labeling of the plot"
    print " -A        Optionally plot the absolute value of the cross correlation"
    print " -D        Optionally print debugging Info"
    print " -O        Optionally swap the I/Q sample time order"
    print " -Z        Optionally Zero event before cross correlation"
    print ""
    print "Glen Langston - 2019 December 10"
    exit()

namearg = 1
iarg = 1          # start searching for input flags

# for all arguments, read list and exit when no flag argument found
while iarg < nargs:
    if sys.argv[iarg].upper() == '-T':   # now look for flags with arguments
        iarg = iarg+1
        mytitle = sys.argv[iarg]
        print 'Plot title: ', mytitle
    elif sys.argv[iarg].upper() == '-D':   # if debugging
        doDebug = True
        print "-D Additional Debug printing."
    elif sys.argv[iarg].upper() == '-Z':   # if zeroing event
        doZero = True
        print "-Z Zeroing event before correlation."
    elif sys.argv[iarg].upper() == '-O':   # if swapping the I/Q time order
        doOrder = True
        print "-O Testing a swap of I/Q sample time order."
    elif sys.argv[iarg].upper() == '-A':   # if absolute value of cross correlation plotted
        doAbs = True
        print "-A Plotting absolute value of cross correlation."
    else:
        break
    iarg = iarg+1
corrarg = iarg  
namearg = corrarg + 1

if namearg >= nargs:
    print "Need a file name to correllate!"
    exit()

print "corrarg: ", sys.argv[corrarg] 
nCorr = int(sys.argv[corrarg])
# keep the cross correllation reasonable
if nCorr < 3:
    nCorr = 3
print "Cross correlating over a range of +/- %d samples" % (nCorr)

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
        print "%s == %s == %12.6f" %( time, timestr, seconds)
    parts = timestr.split('T')
    if len(parts) > 1:
        atime = parts[1]
    else:
        atime = timestr

    if doDebug:
        print "%s == %s" %( time, atime)

    gallon = rs.gallon
    gallat = rs.gallat
    ra = rs.ra
    dec = rs.dec
    label = '%s R,D: %6.2f,%6.2f, Lon,Lat=%5.1f,%5.1f' % ( time,rs.ra,rs.dec,gallon,gallat)
    xs = rs.xdata * 1.E6
    ya = rs.ydataA
    yb = rs.ydataB
    if rs.nTime < 1:
        print "Not an Event: ",filename
        continue

    xv = np.zeros(rs.nSamples*2)
    yv = np.zeros(rs.nSamples*2)
    j = 0
    # convert time offsets to micro-seconds
    dt = 1.E6 * 0.5/rs.bandwidthHz
    t = -2. * dt * rs.refSample
    if iii == -1:  # no op
        print "First Time %12.9f (s); Delta T = %12.9f (s)" % (t, dt)
    for i in range(rs.nSamples):
        if doOrder:
            yv[j] = yb[i]
            yv[j+1] = ya[i]
        else:
            yv[j] = ya[i]
            yv[j+1] = yb[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        xv[j] = t
        j = j + 1
        t = t + dt

    dts = np.zeros( 2*nCorr+1)
    corrs = np.zeros( 2*nCorr+1)

    if doZero:
        center = int(rs.nSamples/2)
        yv[(center-5):(center+5)] = 0.0

    y1s = yv[nCorr:(rs.nSamples-(nCorr+1))]
    if doDebug:
        print "y1s: ", y1s

    printcount = 0
    for kkk in range((2*nCorr)+1):
        # prepare to plot correlation
        jjj = kkk-nCorr
        dts[kkk] = jjj*dt

        y2s = yv[(jjj+nCorr):(jjj + (rs.nSamples-(nCorr+1)))]
        data = pd.DataFrame({'y1': y1s[:], 'y2': y2s[:]})
        zzs = data.corr()
        if doDebug:
            print "zzs: ", zzs
        correlation = zzs.iloc[[0],[1]]
        if doDebug:
            print 'Corr: ', correlation
        corrs[kkk] = np.float(correlation.iloc[0])
        
        if printcount % 50 == 0:
            print "Self correlation offset: %4d, %10.3fus %8.5f" % (jjj, dts[kkk], corrs[kkk])
        printcount += 1

    if doAbs:
        corrs = np.abs( corrs)

    n2 = int( rs.nSamples/2.5)
    ymin = min(corrs)
    ymax = max(corrs)
    yrms = np.std( corrs)

    count = rs.count

    if mytitle == "":
        mytitle = rs.site

    ypeak = max( -ymin, ymax)
    if yrms > 0:
        snr = ypeak/yrms
    else:
        snr = 0.

    print(' Ra: %6.2f Dec: %6.2f Max: %8.3f +/- %7.3f SNR: %6.1f ; %s %s' % (ra, dec, ypeak, yrms, snr, count, label))
    if nplot <= 0:
        title = mytitle + " Az,El: %6.1f,%6.1f" % (rs.telaz, rs.telel)
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(title)
        nplot = nplot + 1
    note = rs.noteA
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.ylim(yallmin,min(yallmax,0.5))


    plt.plot(dts, corrs, colors[iii-1], linestyle=linestyles[iii-1],label=label)
plt.xlim(dts[int(nCorr/2)],dts[2*nCorr])
plt.title(title)
plt.xlabel('Correllation Time Offset (micro-seconds)')
plt.ylabel('Cross Correlation Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
