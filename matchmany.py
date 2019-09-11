#Python Script to plot and the fourier transform of blocks of raw NSF events 
#HISTORY
#19SEP08 GIL match several directories
#19APR15 GIL first working version of event matchin
#19APR14 GIL initial version of event matching
#19APR11 GIL improve labeling
#19APR10 GIL initial version based on FFT
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
iOffset = 1
ifile = 1
iii = ifile
offset = 0.5          # default match offset is  0.5 seconds
nday = 24             # by default divide day in 24 hours
nblock = 256          # number of samples to FFT
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot

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
    print "MATCHMANY: Match events listed  directories"
    print "Usage: MATCH [-OF seconds] dir1 dir2 [.. dirN]"
    print "Where: Optionally the user provides the maximum offset to call a match"
    exit()

    

def main(dirs):
    """
    Main executable for gridding astronomical data
    """
    dpi = 1

    nargs = len(sys.argv)
    if nargs < 2:
        print 'MATCH: MATCH pairs of events in directories'
        print 'usage: MATCH [-OF seconds] dir1 dir2'
        exit()

    dirs = sys.argv[ifile:nargs]
    ndirs = nargs - ifile]

    dir1 = sys.argv[ifile]
    dir2 = sys.argv[ifile+1]
    print "Dir 1: ", dir1
    print "Dir 2: ", dir2
    gridtype = 'PULSAR'
    idir = 0 # index into directory
    from os import listdir
    from os.path import isfile, join

    for dir in dirs:
        files = [f for f in listdir(dir) if isfile(join(dir, f))]
        count = 0
        print "%5d Files in Directory %s" % (len(files), dir)
        events = files
        nEve = 0
        for filename in files:
            parts = filename.split(".")
            nparts = len(parts)
            if nparts < 2:   # if not fooo.eve type file name
                continue
            if parts[nparts-1] == "eve":
                events[nEve] = filename
                nEve = nEve + 1
            else:
                continue

        #end of getting all file names
        #prepare to read values
        mjds = np.zeros(nEve)
        peaks = np.zeros(nEve)
        rmss = np.zeros(nEve)
        dts = np.zeros(nEve)
        iis = np.arange(nEve) * 0
        events = events[0:nEve1]

        for filename in events:
            fullname = join(dir1, filename)
            rs.read_spec_ast( fullname)
            rs.azel2radec()    # compute ra,dec from az,el
            mjds[iii] = rs.emjd
            peaks[iii] = rs.epeak
            rmss[iii] = rs.erms
            if 50 * int(iii/50) == iii:
                print "Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms)
            iii = iii + 1

         print "%5d Events in Directory: %s" % (nEve, dir)

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

    for iii in range(nEve1):
        if dt1s[iii] < offset:
            dts = dt1s[iii]/OneMjdSec
            jjj = ii1s[iii]
            print "Event %s Matches %s; Offset: %9.6f s" % (event1s[iii], event2s[jjj], dts)
            print "%5d %18.9f: %5d %18.9f" % (iii, mjd1s[iii], jjj, mjd2s[jjj])

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

#    for iii in range(nEve2):
#        if dt2s[iii] < offset:
#            dts = dt2s[iii]/OneMjdSec
#            jjj = ii2s[iii]
#            print "Event %s Matches %s; Offset: %9.6f s" % (event2s[iii], event1s[jjj], dts)
#            print "%5d %18.9f: %5d %18.9f" % (iii, mjd2s[iii], jjj, mjd1s[jjj])


if __name__ == "__main__":
    main()

