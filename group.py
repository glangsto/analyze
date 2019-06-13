#Python Script to group events by time
#HISTORY
#19JUN13 GIL initial version
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
eventCounts = np.zeros(nday)  # count events per fraction of a day
eventAveGLat = np.zeros(nday) # count events per fraction of a day
eventAveGLon = np.zeros(nday) # count events per half hour
eventAveRa = np.zeros(nday) # count events per fraction of a day
eventAveDec = np.zeros(nday) # count events per half hour
nu = np.zeros(nblock)  # creat frequency array
nchan = (nblock/2)-1
ysum = np.zeros(nchan)
yfilesum = np.zeros(nchan)
yp2 = np.zeros(nchan)
nsum = 0      # count total number of spectra summed
nFilesum = 0  # count number of spectra in a file
nplot = 0
nFiles = nargs-ifile

print "Number of Files:                ",nFiles
if nFiles < 1:
    print "GROUP: Read all input events and group by specified time interval"
    print "Usage: GROUP [-dt <time interval>] [-sigma <n sigma>] [-note <note for plot title>] <file 1> [<file 2>] ... [<file N>]"
    print "Where optionally the following paramters may be applied"
    print " -dt <time interval> - Time, in seconds, to call a group. Default 60"
    print " -sigma <n sigma>    - Number of sigma of a sample to declare an event"
    print " -note <text>        - Note for the top of the plot"
    exit()
    
print "First File     :                ",sys.argv[ifile]

maxMagnitude = 0.
maxEvent = radioastronomy.Spectrum()
maxFile = ""
#Create a gigantic array of events

events = [ radioastronomy.Spectrum() for i in range(nFiles) ]
eventtimes = np.zeros(nFiles)
eventCount = 0

# now for all files, read and store
for iii in range(nFiles):

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
    
    yv = np.zeros(2*rs.nSamples)
    yc = np.zeros(rs.nSamples,dtype=np.complex)
    j = 0
    for i in range(rs.nSamples):
        yv[j] = rs.ydataA[i]
        yc[i] = rs.ydataA[i] + 1j*rs.ydataB[i]
        j = j + 1
        yv[j] = rs.ydataB[i]
        j = j + 1

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
        ipart = np.int(seconds/(86400./nday))
    #    print "Utc (%8.1f): %s %2d" % (seconds, rs.utc, half)
        eventCounts[ipart] += 1
        eventAveGLat[ipart] += rs.gallat
        eventAveGLon[ipart] += rs.gallon
        eventAveRa[ipart] += rs.ra
        eventAveDec[ipart] += rs.dec

    events[eventCount] = copy.deepcopy( rs)
    eventCount = eventCount + 1

note = rs.noteA

# now compute spectrum of maximum event
# 
print "Event Time         : ", maxEvent.utc
print "Event File         : ", maxFile
print "Maximum Event Sigma: %8.2f" % ( maxMagnitude)
print "Event RA           :   %8.4f d Dec %8.4f d" % (maxEvent.ra, maxEvent.dec)
print "Event G Lon        :   %6.2f d  Lat %6.2f d" % (maxEvent.gallon, maxEvent.gallat)

datetime = "%s" % ( maxEvent.utc)
parts = datetime.split('.')
date = parts[0]

# now summarize events per ipart hour of observation

ntotal = 0
nhours = 0
for ipart in range(nday):
    ntotal += np.float(eventCounts[ipart])
    if eventCounts[ipart] > 0:
        nhours += 1
        
print ""
print "# Total Event  Count:  ", ntotal
print "# Hours Event  G Lon   G Lat    Ra      Dec"
print "# (Utc) Count   (d)     (d)     (h)     (d)"
for ipart in range(nday):
    if eventCounts[ipart] > 0:
        eventAveGLon[ipart] = eventAveGLon[ipart]/np.float(eventCounts[ipart])
        eventAveGLat[ipart] = eventAveGLat[ipart]/np.float(eventCounts[ipart])
        eventAveRa[ipart] = eventAveRa[ipart]/np.float(eventCounts[ipart])
        eventAveDec[ipart] = eventAveDec[ipart]/np.float(eventCounts[ipart])
        print "%6.3f %5d %7.2f %7.2f %7.2f %7.2f" % ( ipart*24./nday, eventCounts[ipart], eventAveGLon[ipart], eventAveGLat[ipart], eventAveRa[ipart]/15., eventAveDec[ipart])
    else:
        if nhours > 5:  # if several hours have data, print zeros to simpiify plotting
            print "%6.3f %5d     0       0       0       0 " % ( ipart*24./nday, eventCounts[ipart])
