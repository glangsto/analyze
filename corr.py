
#Python Script to correlate samples in raw NSF events 
#plot the self correlation values samples near an event
#HISTORY
#21AUG23 GIL add diagnostics
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
doExclude = False
nChannel = 0
printcount = 0

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
    print(" -X        Exclude self cross correllations")
    print(" -Z        Optionally Zero event before cross correlation")
    print("")
    print("Glen Langston - 2019 December 11")
    exit()

nCenter = 0
namearg = 1
iarg = 1          # start searching for input flags

def xcorr( center, nSample, nCorr, yia, yqa, yib, yqb):
    """
    Cross correlation a pair of i,q values
    The input IA and QA are not shifted and are a subset of
    of all samples available.  
    IB and QB are the entire event samples and sub-sets are selected
    inputs:
    center:  Offset from center cross correllation
    nSample: Total (max) number of samples in event)
    nCorr    Half size of cross correlation
    yia, yqa:A event to cross correlate
    yib, yqb:B event to cross correlate
    """
    printcount = 0
    center = int(center)
    nSample = int(nSample)
    n2 = int(nSample/2)
    nCorr = int(nCorr)
    nxcor = (2*nCorr)+1
    xcori = np.zeros(nxcor)
    xcorq = np.zeros(nxcor)
    xsum  = np.zeros(nxcor)

    if doDebug:
        print("center, nSample, nCorr: %d, %d, %d: " % (center, nSample, nCorr))

    yi1 = yia[nCorr:nSample-nCorr-1]
    yq1 = yqa[nCorr:nSample-nCorr-1]
    # for all samples to cross correllate
    for kkk in range(nxcor):
        # select range for cross correllation
        ib = kkk
        ie = kkk + nSample - nxcor

# step through different correlation offsets
        yi2 = yib[ib:ie]
        yq2 = yqb[ib:ie]
        if doMag:
            yi2 = yi2*yi2
            yq2 = yq2*yq2
            sumi = yi2.sum()
            sumq = yq2.sum()
        else:
            yix = yi1*yi2
            yqx = yq1*yq2
            sumi = yix.sum()
            sumq = yqx.sum()
        xcori[kkk] = sumi
        xcorq[kkk] = sumq
        xsum[kkk] = sumi - sumq

        # simultaneously multiple the two time vectors togethere
    if doDebug:
        print( "Size yia, yqa: %d, %d" % (len(yia), len(yqa)))
        print( "Size y1b, yqb: %d, %d" % (len(yib), len(yqb)))
# end of xcorr()
    return xsum

# for all arguments, read list and exit when no flag argument found
while iarg < nargs:
    if sys.argv[iarg].upper() == '-T':   # now look for flags with arguments
        iarg = iarg+1
        mytitle = sys.argv[iarg]
        print('Plot title: ', mytitle)
    elif sys.argv[iarg].upper() == '-X':   # if debugging
        doExclude = True
        print("-X Excluding plots of an antenna with itself.")
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
nread = 0

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
    namea = filename
    nameparts = namea.split("-")
    pia = nameparts[0]
    timea = time
    xs = rs.xdata * 1.E6
    if rs.nTime < 1:
        print("Not an Event: ",filename)
        continue

    center = int(rs.nSamples/2)
    if doZero:
        ns2 = rs.nSamples/2
        print( "Zeroing Center sample: %d+/-5 of %d " % (ns2, rs.nSamples))
        rs.ydataA[(ns2-5):(ns2+5)] = 0.0
        rs.ydataB[(ns2-5):(ns2+5)] = 0.0

    nsamp = rs.nSamples
    # optionally add an offset of samples
    idelta = 0
    # if first plot, must auto-correlate
    yib = rs.ydataA[0:nsamp]
    yqb = rs.ydataB[0:nsamp]
    n = nsamp
    yiarms = np.std( yib)
    yqarms = np.std( yqb)
    print("I rms: %7.4f;  Q rms: %7.4f;" % ( yiarms, yqarms))
    nread = nread + 1
    if nread < 2:
        print( "First file read: n = %d, %d" % (len(yib), len(yqb)))
        # save samples for next correlation
        yifirst = yib
        yqfirst = yqb
        pifirst = pia
        timefirst = timea
        # convert time offsets to micro-seconds
        dt = 1.E6/rs.bandwidthHz
        t = - dt * nCorr
        if iii == -1:  # no op
            print("First Time %12.9f (s); Delta T = %12.9f (s)" % (t, dt))
        # prepare to plot correlation time axis
        dts = np.zeros( (2*nCorr)+1)
        for i in range((2*nCorr)+1):
            dts[i] = (i-(nCorr+1))*dt
        yia = yib
        yqa = yqb
        nameb = namea
        pib = pia
        timeb = timea
        continue
    
    else:
        label = "%s %s X %s %s" % (pia, timea, pib, timeb)
        # maximum correlation comes from the auto correlation
        corrs = xcorr( center, nChannel, nCorr, yia, yqa, yib, yqb)
        ymax = corrs.max()
        ymin = corrs.min()
        yrms1 = np.std( corrs[0:int(nCorr/3)])
        yrms2 = np.std( corrs[int(2*nCorr/3):int(2*nCorr)])
        yrms = (yrms1 + yrms2)/2.
        print("Corr Max: %7.3f, Min: %7.3f Rms: %6.3f " % (ymax, ymin, yrms))
        if mytitle == "":
            mytitle = rs.site
        # prepare for next cross

        ypeak = max( -ymin, ymax)
        if yrms > 0:
            snr = ypeak/yrms
        else:
            snr = 0.

        if nplot <= 0:
            title = mytitle + " Az,El: %6.1f,%6.1f" % (rs.telaz, rs.telel)
            fig,ax1 = plt.subplots(figsize=(10,6))
            fig.canvas.set_window_title(title)
            nplot = 1
        note = rs.noteA
        yallmin = min(ymin,yallmin)
        yallmax = max(ymax,yallmax)
#    plt.ylim(yallmin,yallmax)

        if (not doExclude) or (pia != pib):
            plt.plot(dts, corrs, colors[nplot], \
                linestyle=linestyles[nplot],label=label)
            nplot = nplot + 1

# now prepare for the next correlation    
    nameb = namea
    pib = pia
    timeb = timea
    yia = yib
    yqa = yqb

# now final cross correlation is between first and last
# if more than 2 read
if nread > 2:
        # maximum correlation comes from the auto correlation
    label = "%s %s X %s %s" % (pia, timea, pifirst, timefirst)
    if (not doExclude) or (pia != pib):
        corrs = xcorr( center, nChannel, nCorr, yia, yqa, yifirst, yqfirst)
        plt.plot(dts, corrs, colors[nplot], \
             linestyle=linestyles[nplot],label=label)

plt.xlim(dts[0],dts[2*nCorr])
plt.title(title)
plt.xlabel('Correllation Time Offset (micro-seconds)')
plt.ylabel('Cross Correlation Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
