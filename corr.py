#Python Script to correlate samples in raw NSF events 
#plot the self correlation values samples near an event
#HISTORY
#20FEB09 GIL update for python3
#19DEC13 GIL add magnitude option
#19DEC12 GIL try to signal from cross correlation
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
doZero = False
doOrder = False
nChannel = 0

if nargs < 2:
    print("CORR: Cross correlate Event Samples over a range in sample offsets")
    print("Usage: ")
    print(" CORR [-T 'My Plot Title'] [-A -D -Z] [-C <n>] <nSamples> fileNames ")
    print("where:")
    print(" <nSamples> is the maximum number of samples to compute the cross correlation")
    print("  A range of +/- nSamples correllations are performed and plotted")
    print("  -T <plot title String> Optionally Enables user labeling of the plot")
    print(" -A        Optionally plot the absolute value of the cross correlation")
    print(" -C <n>    Optionally only compute the correlation for the central <n> samples around an event.")
    print(" -D        Optionally print debugging Info")
    print(" -Z        Optionally Zero event before cross correlation")
    print("")
    print("Glen Langston - 2019 December 11")
    exit()

nCenter = 0
namearg = 1
iarg = 1          # start searching for input flags

def xcorr( center, nChannel, nCorr, yi1, yq1, yib, yqb):
    """
    Cross correlation a pair of i,q values
    """
    printcount = 0
    center = int(center)
    nChannel = int(nChannel)
    nCorr = int(nCorr)

    # for all samples to cross correllate
    for kkk in range((2*nCorr)+1):
        jjj = kkk-nCorr
        # select range for cross correllation
        if (center - nChannel + jjj) <= 0:
            ib = 0
            ie = (2*nCorr) + 1
        else:
            ib = center - nChannel + jjj
            ie = center + nChannel + jjj

# step through different correlation offsets
        yi2 = yib[ib:ie]
        yq2 = yqb[ib:ie]
        if doMag:
            yi2 = yi2*yi2
            yq2 = yq2*yq2
        y2 = pd.DataFrame({ 'i': yi2, 'q': yq2})
        if doDebug:
            print("y2.corr(): ", y2.corr())

        # simultaneously multiple the two time vectors togethere
        allcorr = (yi1*yi2) + (yq1*yq2)
        allcorr = allcorr/ymax
        if doDebug:
            print(allcorr)
        corrrms = np.std(allcorr[0:10])
        sumall = allcorr.sum()
        corrs[kkk] = sumall
# end of xcorr()
    return corrs


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
    elif sys.argv[iarg].upper() == '-A':   # if absolute value of cross correlation plotted
        doAbs = True
        print("-A Plotting absolute value of cross correlation.")
    elif sys.argv[iarg].upper() == '-M':   # if correlating magnitudes
        doMag = True
        print("-M Correllating sqrt(I^2 + Q^2) intensity")
    elif sys.argv[iarg].upper() == '-C':   # setting a range of samples to correllate
        iarg = iarg+1
        nChannel = int( sys.argv[iarg])
        print("Summing only over the central +/- %d channels" % (nChannel))
    else:
        break
    iarg = iarg+1
corrarg = iarg  
namearg = corrarg + 1
nnames = nargs - namearg

if nnames <= 1:
    print("Need a file name to correllate!")
    exit()

nnames = min(nnames, 25)  # avoid processing too many files

print("corrarg: ", sys.argv[corrarg]) 
nCorr = int(sys.argv[corrarg])
# keep the cross correllation reasonable
if nCorr < 3:
    nCorr = 3
print("Cross correlating over a range of +/- %d samples" % (nCorr))

note = "Events"
nplot = 0

iname = namearg   # prepare to read the first file

rs = radioastronomy.Spectrum()

for iii in range(nnames):
    filename = sys.argv[iname]
    iname = iname + 1

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

    center = int(rs.nSamples/2)
    if doZero:
        print( "Zeroing Center sample: %d+/-5 %d " % (center, rs.nSamples))
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
    # optionally add an offset of samples
    idelta = 9
    # if first plot, must auto-correlate
    if nplot == 0:
        yib = rs.ydataA[(ib+idelta):(ie+idelta)]
        yqb = rs.ydataB[(ib+idelta):(ie+idelta)]
        yifirst = yib
        yqfirst = yqb
        yi1 = yib
        yq1 = yqb
    else: # use last event to cross correlate
        yi1 = yib
        yq1 = yqb
        # now save current file for next cross correllate
        yib = rs.ydataA[(ib+idelta):(ie+idelta)]
        yqb = rs.ydataB[(ib+idelta):(ie+idelta)]

    yi1rms = np.std( yi1)
    yq1rms = np.std( yq1)

    y1 = pd.DataFrame({ 'i': yi1, 'q': yq1})
    print("I rms: %7.4f;  Q rms: %7.4f; n: %d" % ( yi1rms, yq1rms, n))
    #    print "y1.corr(): ", y1.corr()

    # convert time offsets to micro-seconds
    dt = 1.E6/rs.bandwidthHz
    t = - dt * rs.refSample
    if iii == -1:  # no op
        print("First Time %12.9f (s); Delta T = %12.9f (s)" % (t, dt))
    for i in range(n):
        xv[i] = t
        t = t + dt
        # prepare to plot correlation time axis

    dts = np.zeros( (2*nCorr)+1)
    for i in range((2*nCorr)+1):
        dts[i] = (i-(nCorr+1))*dt

    # maximum correlation comes from the auto correlation
    allcorr = (yi1*yi1) + (yq1*yq1)
    ymax = allcorr.sum()
    if doMag:
        yi1 = yi1*yi1
        yq1 = yq1*yq1

    corrs = xcorr( center, nChannel, nCorr, yi1, yq1, yib, yqb)

    # all cross correlations complete
    ymax = max(corrs)
    # if plotting absolute values of correllations
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
    plt.ylim(yallmin,yallmax)

    if (printcount % 10 == 0) or ((jjj > -3) and (jjj < 3)):
        print("Self correlation offset: %4d, %10.3fus %8.5f +/- %7.5f" % (jjj, dts[kkk], corrs[kkk], corrrms))
        printcount += 1


    plt.plot(dts, corrs, colors[iii], linestyle=linestyles[iii],label=label)
plt.xlim(dts[0],dts[2*nCorr])
plt.title(title)
plt.xlabel('Correllation Time Offset (micro-seconds)')
plt.ylabel('Cross Correlation Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
