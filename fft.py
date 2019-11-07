#Python Script to plot the fourier transform of blocks of raw NSF events 
#plot the FFT of raw data from an event
#HISTORY
#19FEB27 GIL reduce printout if only a few hour ranges have data
#19FEB25 GIL compute average locations
#19FEB21 GIL initial version
#
import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
from scipy.fftpack import fft
import numpy as np
import copy
from scipy.signal import blackman

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
ncolor = min( len(colors), 7)

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

# separate arguments form file names
ifile = 1
iii = ifile

nday = 24             # by default divide day in 24 hours
nblock = 256          # number of samples to FFT
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot
# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-BL":
        nblock = np.int( sys.argv[iii+1])
        iii = iii + 1
        print "FFT Block Size: ", nblock
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-KP":
        kpercount = np.float( sys.argv[iii+1])
        iii = iii + 1
        print "Kelvins Per Count: ", kpercount
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-ND":
        nday = np.int( sys.argv[iii+1])
        iii = iii + 1
        print "Divide Day into N Parts: ", nday
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-SI":
        sigma = np.float( sys.argv[iii+1])
        iii = iii + 1
        print "Keeping Events > %7.2f Sigma " % (sigma)
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print "Note: ", note
        ifile = ifile + 2
        aFix = True
    iii = iii + 1

N = nblock                 # abreviation
ablock = np.zeros(nblock)  # creat array to FFT
eventblock = np.zeros(nblock)  # creat array to FFT
eventCounts = np.zeros(nday)  # count events per fraction of a day
eventAveGLat = np.zeros(nday) # count events per fraction of a day
eventAveGLon = np.zeros(nday) # count events per half hour
eventAveRa = np.zeros(nday) # count events per fraction of a day
eventAveDec = np.zeros(nday) # count events per half hour
w = blackman(N)
nu = np.zeros(nblock)  # creat frequency array
nchan = (nblock/2)-1
ysum = np.zeros(nchan)
yfilesum = np.zeros(nchan)
yp2 = np.zeros(nchan)
nsum = 0      # count total number of spectra summed
nfilesum = 0  # count number of spectra in a file
nplot = 0
nfiles = nargs-ifile

print "Number of Files:                ",nfiles
if nfiles < 1:
    print "FFT: Fourier Transform and sum a time series of events"
    print "Usage: FFT [-bl <n samples>] [-sigma <n sigma>] [-nd <n day parts>] [-kp <kelvins per count>] [-note <note for plot title>] <file 1> [<file 2>] ... [<file N>]"
    print "Where optionally the following paramters may be applied"
    print " -bl <n samples>     - Number of samples to FFT to produce a spectra (power of 2"
    print " -sigma <n sigma>    - Number of sigma of a sample to declare an event"
    print " -kpercount <factor> - Gain factor to convert counts to Kelvin"
    print "   Estimate by an assumed system temperature for the band pass"
    print "   And running FFT without this factor to get the value in counts"
    print " -note <text>        - Note for the top of the plot"
    exit()
    
print "First File     :                ",sys.argv[ifile]

maxMagnitude = 0.
maxEvent = radioastronomy.Spectrum()
maxFile = ""

for iii in range(nfiles):

    filename = sys.argv[iii+ifile]

    rs = radioastronomy.Spectrum()
#    print filename
    rs.read_spec_ast( filename)
    if note != "":
        rs.noteA = note
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
    xs = rs.xdata 
    if rs.nTime < 1:
        print "Not an Event: ",filename
        continue

    xv = np.zeros(rs.nSamples*2)
    yv = np.zeros(rs.nSamples*2)
    yc = np.zeros(rs.nSamples,dtype=np.complex)
    ymag = np.zeros(rs.nSamples)
    nSamples = 2L * rs.nSamples
    j = 0
    dt = 0.5/rs.bandwidthHz
    t = xs[0]
    for i in range(rs.nSamples):
        yv[j] = rs.ydataA[i]
        yc[i] = rs.ydataA[i] + 1j*rs.ydataB[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        yv[j] = rs.ydataB[i]
        xv[j] = t
        j = j + 1
        t = t + dt

    # compute magnitude of I/Q samples with vector math
    ymag = np.absolute(yc)
    # now find the maximum event in data series
    ymagmax = max(ymag)
    # compute RMS of entire time series
    yrms = np.sqrt(yv.dot(yv)/yv.size)
    # if a valid, noisy observation
    if yrms > 0. :
        nsigma = ymagmax/yrms
        if nsigma > maxMagnitude:
            maxEvent = copy.deepcopy(rs)
            maxMagnitude = nsigma
            maxFile = filename
    else:
        print "Problem observation: %7.3f rms in file %s\n", yrms, filename
        continue
    # if not a signficant event

    if nsigma > sigma:
    # now compute time averages of number of events and locations
        midnight = rs.utc.replace(hour=0, minute=0, second=0, microsecond=0)
        seconds = (rs.utc - midnight).seconds
        iday = np.int(seconds/(86400./nday))
    #    print "Utc (%8.1f): %s %2d" % (seconds, rs.utc, half)
        eventCounts[iday] += 1
        eventAveGLat[iday] += rs.gallat
        eventAveGLon[iday] += rs.gallon
        eventAveRa[iday] += rs.ra
        eventAveDec[iday] += rs.dec

    # will compute several FFTs, discarding any extra samples at end
    nfft = nSamples/nblock
    lll = 0
    # 
    BW = 1.E-6*rs.bandwidthHz
    # frequency axis is always the same
    nu = np.linspace(0.0, BW, nblock/2) + (1.E-6*(rs.centerFreqHz-rs.bandwidthHz/2.))
    xmin = min(nu)
    xmax = max(nu)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)
    # zero sum of all spectra in this file

    yfilesum = np.zeros(nchan)
    nfilesum = 0  # count number of spectra in a file

    for jjj in range(nfft):
        # for each block of samples
        for kkk in range(nblock):
            ablock[kkk] = yv[lll]
            lll += 1
        yf = fft(ablock*w)  # apply blackman window/weighting to samples
#        yf = fft(ablock)   # fft without window
        xp = nu[1:N//2]
        yp = 2.0/N*np.abs(yf[1:N//2])
        ysum = ysum + yp
        nsum += 1
        # average all spectra in this file.
        yfilesum = yfilesum + yp
        nfilesum += 1

        ymin = min(yp)
        ymax = max(yp)

        if jjj == 0:
            print(' Max: %9.3f Min: %9.3f ; %3d %s' % (ymax*kpercount, ymax*kpercount, nfft, label))
#        plabel = ("%d" % (jjj))
    if nplot <= 0:
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(date)
        nplot = nplot + 1

    # prepare to plot sum of all spectra in this file
    if nfilesum > 0:
        yp = yfilesum/np.float(nfilesum)
    else:
        print "problem with sum of spectra in a afile: ",filename
        continue
    # flip spectra
    nnn = (nchan/2)
    for kkk in range(nchan):
        yp2[kkk] = yp[nnn]
        nnn -= 1
        if nnn < 0:
            nnn = nchan - 1
            
#    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)

# plot the middle spectrum for each events, until out of colors
    if jjj != nfft-1:
        continue
    if nplot >= ncolor:
        continue
    nplot = nplot + 1
    yp2 = kpercount * yp2
    label = '%s %s Lon,Lat=%5.1f,%5.1f' % ( date, time, gallon, gallat)
    plt.plot(xp, yp2, colors[nplot-1], linestyle=linestyles[nplot-1],label=label)
note = rs.noteA
print "Number of FFTs per time series: ",nfft
print "Number of spectra summed:       ",nsum
print "Seconds of observations:        ",nsum*nblock*1.e-6/(2.*BW)

# flip spectra
ysum = ysum / np.float(nsum)
nnn = (nchan/2)
for kkk in range(nchan):
    yp2[kkk] = ysum[nnn]
    nnn -= 1
    if nnn < 0:
        nnn = nchan - 1

yp2 = kpercount * yp2

plabel = "Sum of %5d" % nsum
# finally plot the sum of all spectra; show band pass accuracy 
plt.plot(xp, yp2, colors[0], linestyle=linestyles[0],label=plabel, lw=2)
#print 'XP: ', xp
#print 'YP: ', yp

# now compute spectrum of maximum event
# 
print "Event Time         : ", maxEvent.utc
print "Event File         : ", maxFile
print "Maximum Event Sigma: %8.2f" % ( maxMagnitude)
print "Event RA           :   %8.4f d Dec %8.4f d" % (maxEvent.ra, maxEvent.dec)
print "Event G Lon        :   %6.2f d  Lat %6.2f d" % (maxEvent.gallon, maxEvent.gallat)

# event is in the middle of the data samples
xv = np.zeros(maxEvent.nSamples*2)
yv = np.zeros(maxEvent.nSamples*2)
ymag = np.zeros(maxEvent.nSamples)
nSamples = 2L * maxEvent.nSamples
j = 0
dt = 0.5/maxEvent.bandwidthHz
t = xs[0]
for i in range(maxEvent.nSamples):
    yv[j] = maxEvent.ydataA[i]
    xv[j] = t
    j = j + 1
    t = t + dt
    yv[j] = maxEvent.ydataB[i]
    xv[j] = t
    j = j + 1
    t = t + dt

# center the event in the FFT window    
lll = (maxEvent.nSamples + nblock)/2

# for each block of samples
for kkk in range(nblock):
    ablock[kkk] = yv[lll]
    lll += 1
yf = fft(ablock)
xp = nu[1:N//2]
yp = 2.0/N*np.abs(yf[1:N//2])

# flip spectra
nnn = (nchan/2)
for kkk in range(nchan):
    yp2[kkk] = yp[nnn]
    nnn -= 1
    if nnn < 0:
        nnn = nchan - 1

yp2 = kpercount * yp2

datetime = "%s" % ( maxEvent.utc)
parts = datetime.split('.')
date = parts[0]

# finally plot the FFT of the event
plt.plot(xp, yp2, colors[8], linestyle=linestyles[8],label="Event "+date, lw=3)

plt.xlim(xallmin,xallmax)
plt.title(note)
plt.xlabel('Frequency (MHz)')
if kpercount == 1.:
    plt.ylabel('Intensity (Counts)')
    plt.ylim(yallmin,1.25*yallmax)
else:
    plt.ylabel('Intensity (Kelvins)')
    plt.ylim(kpercount*yallmin,kpercount*1.25*yallmax)
plt.legend(loc='upper right')
plt.show()

# now summarize events per iday hour of observation

ntotal = 0
nhours = 0
for iday in range(nday):
    ntotal += np.float(eventCounts[iday])
    if eventCounts[iday] > 0:
        nhours += 1
        
print ""
print "# Total Event  Count:  ", ntotal
print "# Hours Event  G Lon   G Lat    Ra      Dec"
print "# (Utc) Count   (d)     (d)     (h)     (d)"
for iday in range(nday):
    if eventCounts[iday] > 0:
        eventAveGLon[iday] = eventAveGLon[iday]/np.float(eventCounts[iday])
        eventAveGLat[iday] = eventAveGLat[iday]/np.float(eventCounts[iday])
        eventAveRa[iday] = eventAveRa[iday]/np.float(eventCounts[iday])
        eventAveDec[iday] = eventAveDec[iday]/np.float(eventCounts[iday])
        print "%6.3f %5d %7.2f %7.2f %7.2f %7.2f" % ( iday*24./nday, eventCounts[iday], eventAveGLon[iday], eventAveGLat[iday], eventAveRa[iday]/15., eventAveDec[iday])
    else:
        if nhours > 5:  # if several hours have data, print zeros to simpiify plotting
            print "%6.3f %5d     0       0       0       0 " % ( iday*24./nday, eventCounts[iday])
