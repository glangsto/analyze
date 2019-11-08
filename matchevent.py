#Python Script to plot and the fourier transform of blocks of raw NSF events 
#HISTORY
#19NOV08 GIL reduce printout
#19OCT10 GIL optionall plot matchs
#19APR15 GIL first working version of event matchin
#19APR14 GIL initial version of event matching
#19APR11 GIL improve labeling
#19APR10 GIL initial version based on FFT
#
import matplotlib.pyplot as plt
import sys
import os
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
iOffset = 1
ifile = 1
iii = ifile
offset = 0.5          # default match offset is  0.5 seconds
nday = 24             # by default divide day in 24 hours
nblock = 256          # number of samples to FFT
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot
doPlot = False
doDebug = False

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-OF":
        offset = np.float( sys.argv[iii+1])
        iii = iii + 1
        print "Maximum Time Offset: %8.6f s for Match: " % ( offset)
        offset = offset/86400.   # convert to MJDs
        ifile = ifile + 2
    if str(anarg[0:3]) == "-ND":
        ndays = np.int( sys.argv[iii+1])
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
    if anarg[0:2] == "-P":
        doPlot = True
        print "Plotting Matching events"
        ifile = ifile + 1
    if anarg[0:2] == "-D":
        doDebug = True
        print "Debugging"
        ifile = ifile + 1
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
nplot = 0
nfiles = nargs-ifile

MAXEVENTS = 5000

if nfiles < 2:
    print "MATCH: Match events listed in Two directories"
    print "Usage: MATCH [-OF seconds] [-D] dir1 dir2"
    print "Where: Optionally the user provides the maximum offset to call a match"
    print "Where -D  Optionally print debugging info"
    print ""
    print "Glen Langston, November 8, 2019"
    exit()
    

def main():
    """
    Main executable for gridding astronomical data
    """
    dpi = 1

    nargs = len(sys.argv)
    if nargs < 2:
        print 'MATCH: MATCH pairs of events in directories'
        print 'usage: MATCH [-OF seconds] dir1 dir2'
        exit()

    dir1 = sys.argv[ifile]
    dir2 = sys.argv[ifile+1]
    if doDebug:
        print "Dir 1: ", dir1
        print "Dir 2: ", dir2
    gridtype = 'PULSAR'

    from os import listdir
    from os.path import isfile, join
    files1 = [f for f in listdir(dir1) if isfile(join(dir1, f))]
    files2 = [f for f in listdir(dir2) if isfile(join(dir2, f))]

    count = 0
    if doDebug: 
        print "%5d Files in Directory 1" % (len(files1))
        print "%5d Files in Directory 2" % (len(files2))
    nEve1 = 0
    nEve2 = 0

    rs = radioastronomy.Spectrum()

    event1s = files1
    event2s = files2

    for filename in files1:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "eve":
            event1s[nEve1] = filename
            nEve1 = nEve1 + 1
        else:
            continue

    mjd1s = np.zeros(nEve1)
    peak1s = np.zeros(nEve1)
    rms1s = np.zeros(nEve1)
    dt1s = np.zeros(nEve1)
    ii1s = np.arange(nEve1) * 0
    event1s = event1s[0:nEve1]

    iii = 0
    for filename in event1s:
        fullname = join(dir1, filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        mjd1s[iii] = rs.emjd
        peak1s[iii] = rs.epeak
        rms1s[iii] = rs.erms
        if (50 * int(iii/50) == iii) and doDebug:
            print "Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms)
        iii = iii + 1

    for filename in files2:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "eve":
            event2s[nEve2] = filename
            nEve2 = nEve2 + 1
        else:
            continue

    mjd2s = np.zeros(nEve2)
    peak2s = np.zeros(nEve2)
    rms2s = np.zeros(nEve2)
    dt2s = np.zeros(nEve2)
    ii2s = np.arange(nEve2) * 0
    event2s = event2s[0:nEve2]

    iii = 0
    for filename in event2s:
        fullname = join(dir2, filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        mjd2s[iii] = rs.emjd
        peak2s[iii] = rs.epeak
        rms2s[iii] = rs.erms
        if (50 * int(iii/50) == iii) and doDebug:
            print "Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms)
        iii = iii + 1

    if doDebug:
        print "%5d Events in Directory: %s" % (nEve1, dir1)
        print "%5d Events in Directory: %s" % (nEve2, dir2)

    # now match event times
    for iii in range(nEve1):
        mjd1 = mjd1s[iii]
        dtj = 0
        dtMin = abs(mjd1 - mjd2s[0])
        for jjj in range(nEve2):
            mjd2 = mjd2s[jjj]
            dt = abs(mjd1 - mjd2)
            if dt < dtMin:
                dtMin = dt
                dtj = jjj
        ii1s[iii] = dtj
        dt1s[iii] = dtMin

    # now report out the minimum time offsets and times less than 1 second off
    OneMjdSec = 1./86400.
    TwoMjdSec = 2.*OneMjdSec

    nMatch = 0
    for iii in range(nEve1):
        if dt1s[iii] < offset:
            nMatch = nMatch + 1
            dts = dt1s[iii]/OneMjdSec
            jjj = ii1s[iii]
            print "Event %s Matches %s; Offset: %9.6f s" % (event1s[iii], event2s[jjj], dts)
            if doDebug:
                print "%5d %18.9f: %5d %18.9f" % (iii, mjd1s[iii], jjj, mjd2s[jjj])
            if doPlot:
                eventAName = dir1 + "/" + event1s[iii]
                eventBName = dir2 + "/" + event2s[jjj]
                plotEvent = "~/bin/E %s %s" % (eventAName, eventBName)
                os.system(plotEvent)

    print "Found %d Event Matches" % (nMatch)
    # now match event times
    for iii in range(nEve2):
        mjd2 = mjd2s[iii]
        dtj = 0
        dtMin = abs(mjd2 - mjd1s[0])
        for jjj in range(nEve1):
            mjd1 = mjd1s[jjj]
            dt = abs(mjd1 - mjd2)
            if dt < dtMin:
                dtMin = dt
                dtj = jjj
        ii2s[iii] = dtj
        dt2s[iii] = dtMin

if __name__ == "__main__":
    main()

