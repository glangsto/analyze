#Python Script to correlate samples in raw NSF events 
#plot the cross correlation values samples from two horns
#HISTORY
#19DEC13 GIL add magnitude option
#19DEC12 GIL try to estimate signal from cross correlation
#19DEC11 GIL simultaneously correlate I and Q parts 
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
doMag = False
doOrder = False
nChannel = 0
doZero = False

linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

if nargs < 2:
    print("CORR: Cross correlate Event Samples over a range in sample offsets")
    print("Usage: ")
    print(" CORR [-T 'My Plot Title'] [-A -D -M -Z] [-C <n>] <nSamples> fileNames ")
    print("where:")
    print(" <nSamples> is the maximum number of samples to compute the cross correlation")
    print("  A range of +/- nSamples correllations are performed and plotted")
    print("  -T <plot title String> Optionally Enables user labeling of the plot")
    print(" -A        Optionally plot the absolute value of the cross correlation")
    print(" -C <n>    Optionally only compute the correlation for the central <n> samples around an event.")
    print(" -D        Optionally print debugging Info")
    print(" -M        Optionally correlate the Magnitude of the signals")
    print(" -Z        Optionally Zero event before cross correlation")
    print("")
    print("Glen Langston - 2019 December 11")
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
    elif sys.argv[iarg].upper() == '-Z':   # if zeroing event
        doZero = True
        print("-Z Zeroing event before correlation.")
    elif sys.argv[iarg].upper() == '-M':   # if correlating magnitudes
        doMag = True
        print("-M Correllating sqrt(I^2 + Q^2) intensity")
    elif sys.argv[iarg].upper() == '-C':   # setting a range of samples to correllate
        iarg = iarg+1
        nChannel = int( sys.argv[iarg])
        print("Summing only over the central +/- %d channels" % (nChannel))
    elif sys.argv[iarg].upper() == '-A':   # if absolute value of cross correlation plotted
        doAbs = True
        print("-A Plotting absolute value of cross correlation.")
    else:
        break
    iarg = iarg+1
corrarg = iarg  
namearg = corrarg + 1

if namearg >= nargs:
    print("Need a file name to correllate!")
    exit()

print("corrarg: ", sys.argv[corrarg]) 
nCorr = int(sys.argv[corrarg])
# keep the cross correllation reasonable
if nCorr < 3:
    nCorr = 3
print("Cross correlating over a range of +/- %d samples" % (nCorr))

note = "Events"
nplot = 0
count = 0
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
    label = '%s R,D: %6.2f,%6.2f, Lon,Lat=%5.1f,%5.1f' % ( time,rs.ra,rs.dec,gallon,gallat)
    xs = rs.xdata * 1.E6
    if rs.nTime < 1:
        print("Not an Event: ",filename)
        continue

    # count number of event read
    count = count + 1
    center = int(rs.nSamples/2)
    if doZero:
        center = int(rs.nSamples/2)
        rs.ydataA[(center-5):(center+5)] = 0.0
        rs.ydataB[(center-5):(center+5)] = 0.0

    if nChannel <= 0:
        ib = nCorr
        ie = rs.nSamples-nCorr
    else:
        ib = center - nChannel
        ie = center + nChannel
    n  = ie-ib

    if doDebug: 
        print("ib: %d, ie: %d, n: %d" % ( ib, ie, n))
    xv = np.zeros(n)
#    If only have one event can not yet cross correlate
    if count <= 1: 
        yi1 = rs.ydataA[ib:ie]
        yq1 = rs.ydataB[ib:ie]
        if doMag:
            yi1 = yi1*yi1
            yq1 = yq1*yq1
        continue

    # convert time offsets to micro-seconds
    dt = 1.E6/rs.bandwidthHz
    t = - dt * rs.refSample
    if iii == -1:  # no op
        print("First Time %12.9f (s); Delta T = %12.9f (s)" % (t, dt))
    for i in range(n):
        xv[i] = t
        t = t + dt

    dts = np.zeros( 2*nCorr+1)
    corrs = np.zeros( 2*nCorr+1)

    # maximum correlation comes from the auto correlation
    ymax = 0.

    printcount = 0
    for kkk in range((2*nCorr)+1):
        # prepare to plot correlation
        jjj = kkk-nCorr
        if nChannel <= 0:
            ib = kkk
            ie = kkk+n
        else:
            ib = center - nChannel + jjj
            ie = center + nChannel + jjj
        # prepare to plot correlation time axis
        dts[kkk] = jjj*dt

        yi2 = rs.ydataA[ib:ie]
        yq2 = rs.ydataB[ib:ie]
        if doMag:
            yi2 = yi2*yi2
            yq2 = yq2*yq2
        if doDebug:
            y2 = pd.DataFrame({ 'i': yi2, 'q': yq2})
            print("y2.corr(): ", y2.corr())

        allcorr = (yi1*yi2) + (yq1*yq2)
        corrrms = np.std( allcorr[0:10])
        sumall = allcorr.sum()
        corrs[kkk] = sumall
        
        if (printcount % 10 == 0) or ((jjj > -3) and (jjj < 4)):
            print("Self correlation offset: %4d, %10.3fus %8.5f +/- %7.5f" % (jjj, dts[kkk], corrs[kkk], corrrms))
        printcount += 1
    # if a magnitude correllation: then there is a zero offset
    if doMag:
        ymed = np.median(corrs)
        corrs = corrs - ymed
        
    ymax = max(corrs)
    corrs = corrs/ymax
    if doAbs:
        corrs = np.abs( corrs)

    ymin = min(corrs)

    yrms1 = np.std( corrs[0:int(nCorr/3)])
    yrms2 = np.std( corrs[int(2*nCorr/3):int(2*nCorr)])
    yrms = (yrms1 + yrms2)/2.

    print("Corr Max: %7.3f, Min: %7.3f Rms: %6.3f " % (ymax, ymin, yrms))
    if mytitle == "":
        mytitle = rs.site

    ypeak = max( -ymin, ymax)
    if yrms > 0:
        snr = ypeak/yrms
    else:
        snr = 0.

#   now save most recently read event for cross correlation the next round
    if nChannel <= 0:
        ib = nCorr
        ie = rs.nSamples-nCorr
    else:
        ib = center - nChannel
        ie = center + nChannel

    # transfer to first spectrum for next correllatino 
    yi1 = rs.ydataA[ib:ie]
    yq1 = rs.ydataB[ib:ie]
    if doMag:
        yi1 = yi1*yi1
        yq1 = yq1*yq1

#    print(' Ra: %6.2f Dec: %6.2f Max: %8.3f +/- %7.3f SNR: %6.1f ; %s' % (ra, dec, ypeak, yrms, snr, label))
    if nplot <= 0:
        title = mytitle + " Az,El: %6.1f,%6.1f" % (rs.telaz, rs.telel)
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(title)
        nplot = nplot + 1
    note = rs.noteA
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
#    plt.ylim(yallmin,min(yallmax,0.5))
#    plt.ylim(yallmin,yallmax)

    plt.plot(dts, corrs, colors[iii-1], linestyle=linestyles[iii-1],label=label)
#plt.xlim(dts[0],dts[2*nCorr])
plt.title(title)
plt.xlabel('Correllation Time Offset (micro-seconds)')
plt.ylabel('Cross Correlation Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
