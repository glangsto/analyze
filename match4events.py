#Python find matchs in 4 data directories
#HISTORY
#20FEB21 GIL search odroid directory too
#20FEB08 GIL python 2 to 3; search for different names; allow only 3 directories; add plotting
#20JAN14 GIL fix index boundary issues
#19DEC27 GIL if date is after December 20, use different directory
#19NOV30 GIL fix reading first directory
#19NOV29 GIL set print for number of matches
#19NOV27 GIL now check for multiple matches
#19NOV19 GIL add match of all 4 directories
#19NOV17 GIL get matching of pairs of directories working
#19NOV15 GIL read date 
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
if nargs < 2:
    print("MATCH4: Match events listed in four directories")
    print("Usage: MATCH4 [-OF seconds] [-D] [-C Date] [dir1 dir2 dir3 dir4]")
    print("Where: Optionally the user provides the maximum time offset (secs) to call a match")
    print("Where -D  Optionally print debugging info")
    print("Where -C  Optionally provide a calendar date (ie 19Nov17) instead of directories")
    print("Where -N <n> Optionally print matches when number is equal or greater to <n>")
    print("Where -P Plot matches")
    print("")
    print("Glen Langston, November 18, 2019")
    exit()


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
offset = 1.0/86400.   # default match offset is  1 seconds = 1/86400 of a day
nday = 24             # by default divide day in 24 hours
nblock = 256          # number of samples to FFT
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot
calendar = ""
doPlot = False
doDebug = False
nPrint = 4

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-OF":
        offset = np.float( sys.argv[iii+1])
        iii = iii + 1
        print("Maximum Time Offset: %8.6f s for Match: " % ( offset))
        offset = offset/86400.   # convert to MJDs
        ifile = ifile + 2
    if str(anarg[0:3]) == "-N":
        nPrint = np.int( sys.argv[iii+1])
        iii = iii + 1
        if nPrint <= 0:
            nPrint = 2
        print("Print if %d or more matches" % (nPrint))
        ifile = ifile + 2
    if str(anarg[0:3]) == "-SI":
        sigma = np.float( sys.argv[iii+1])
        iii = iii + 1
        print("Keeping Events > %7.2f Sigma " % (sigma))
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print("Note: ", note)
        ifile = ifile + 2
        aFix = True
    if anarg[0:2] == "-C":
        calendar = sys.argv[iii+1]
        iii = iii + 1
        print("Matching Date: ", calendar)
        ifile = ifile + 2
        aFix = True
    if anarg[0:2] == "-P":
        doPlot = True
        print("Plotting Matching events")
        ifile = ifile + 1
    if anarg[0:2] == "-D":
        doDebug = True
        print("Debugging")
        ifile = ifile + 1
    iii = iii + 1

N = nblock                 # abreviation
ablock = np.zeros(nblock)  # creat array to FFT
eventblock = np.zeros(nblock)  # creat array to FFT
w = blackman(N)
nu = np.zeros(nblock)  # creat frequency array
nplot = 0
nfiles = nargs-ifile
dirs = [ "", "", "", "", "", "", "", "", "", "","", "", "", "", "", "",""]

# now find calendars
if calendar == "":
    dirs[0] = sys.argv[ifile]
    dirs[1] = sys.argv[ifile+1]
    ndir = 2
    if nfiles > 2:
        dir[ndir] = sys.argv[ifile+2]
        ndir+1
    if nfiles > 3:
        dir[ndir] = sys.argv[ifile+3]
        ndir+1
else:
    ndir = 0
    for idir in range(1,10):
        adir = "pi%d-events-" % idir
        adir = adir + calendar
        if os.path.isdir(adir):
            dirs[ndir] = adir
            ndir = ndir + 1
            print("Found directory: %s" % adir)

    for idir in range(1,10):
        adir = "odroid%d-events-" % idir
        adir = adir + calendar
#        print("TestDir: %s" % (adir))
        if os.path.isdir(adir):
            dirs[ndir] = adir
            ndir = ndir + 1
            print("Found directory: %s" % adir)

print( "Found %d directories" % (ndir) )

if doDebug:
    print("Dir 1: ", dirs[0])
    print("Dir 2: ", dirs[1])
    print("Dir 3: ", dirs[2])
    print("Dir 4: ", dirs[3])

MAXEVENTS = 10000

def main():
    """
    Main executable for matching transient events
    """
    dpi = 1

    nargs = len(sys.argv)
    if nargs < 2:
        print('MATCH: MATCH pairs of events in directories')
        print('usage: MATCH [-OF seconds] dir1 dir2')
        exit()

    from os import listdir
    from os.path import isfile, join
    files1 = [f for f in listdir(dirs[0]) if isfile(join(dirs[0], f))]
    files2 = [f for f in listdir(dirs[1]) if isfile(join(dirs[1], f))]
    if len(dirs[2]) > 0:
        files3 = [f for f in listdir(dirs[2]) if isfile(join(dirs[2], f))]
    else:
        files3 = ""
    if len(dirs[3]) > 0: 
        files4 = [f for f in listdir(dirs[3]) if isfile(join(dirs[3], f))]
    else:
        files4 = ""

    count = 0
    if doDebug: 
        print("%5d Files in Directory 1" % (len(files1)))
        print("%5d Files in Directory 2" % (len(files2)))
        print("%5d Files in Directory 3" % (len(files3)))
        print("%5d Files in Directory 4" % (len(files4)))
    nEve1 = 0
    nEve2 = 0
    nEve3 = 0
    nEve4 = 0

    rs = radioastronomy.Spectrum()

# now start reading all the files
    event1s = files1
    event2s = files2
    event3s = files3
    event4s = files4

# first directory
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

    if nEve1 < 1:
        nEve1 = 1
    mjd1s = np.zeros(nEve1)
    peak1s = np.zeros(nEve1)
    rms1s = np.zeros(nEve1)
    dt12s = np.zeros(nEve1)
    ii12s = np.arange(nEve1) * 0
    event1s = event1s[0:nEve1]
    if nEve1 == 1:
        nEve1 = 0

    iii = 0
    for filename in event1s:
        fullname = join(dirs[0], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        mjd1s[iii] = rs.emjd
        peak1s[iii] = rs.epeak
        rms1s[iii] = rs.erms
        if (100 * int(iii/100) == iii) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))
        iii = iii + 1

    print("%5d Events in Directory: %s" % (nEve1, dirs[0]))

# second directory
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

# directory 2
    if nEve2 < 1:
        nEve2 = 1
    mjd2s = np.zeros(nEve2)
    peak2s = np.zeros(nEve2)
    rms2s = np.zeros(nEve2)
    event2s = event2s[0:nEve2]

# matches between 2 and 3
    dt23s = np.zeros(nEve2)
    ii23s = np.arange(nEve2) * 0
# matches between 2 and 4
    dt24s = np.zeros(nEve2)
    ii24s = np.arange(nEve2) * 0

    if nEve2 == 1:
        nEve2 = 0

    iii = 0
    for filename in event2s:
        fullname = join(dirs[1], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        mjd2s[iii] = rs.emjd
        peak2s[iii] = rs.epeak
        rms2s[iii] = rs.erms
        if (100 * int(iii/100) == iii) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))
        iii = iii + 1

    print("%5d Events in Directory: %s" % (nEve2, dirs[1]))

# third directory
    for filename in files3:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "eve":
            event3s[nEve3] = filename
            nEve3 = nEve3 + 1
        else:
            continue

# directory 3
    if nEve3 < 1:
        nEve3 = 1
    mjd3s = np.zeros(nEve3)
    peak3s = np.zeros(nEve3)
    rms3s = np.zeros(nEve3)
    event3s = event3s[0:nEve3]

# matches between 1 and 3
    dt13s = np.zeros(nEve1)
    ii13s = np.arange(nEve1) * 0
# match betweeen 3 and 4
    dt34s = np.zeros(max(nEve3,nEve4))
    ii34s = np.arange(max(nEve3,nEve4)) * 0

    if nEve3 == 1:
        nEve3 = 0

    iii = 0
    for filename in event3s:
        fullname = join(dirs[2], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        mjd3s[iii] = rs.emjd
        peak3s[iii] = rs.epeak
        rms3s[iii] = rs.erms
        if (100 * int(iii/100) == iii) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))
        iii = iii + 1

    print("%5d Events in Directory: %s" % (nEve3, dirs[2]))

# forth directory
    for filename in files4:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "eve":
            event4s[nEve4] = filename
            nEve4 = nEve4 + 1
        else:
            continue

# directory 4
    if nEve4 < 1:
        nEve4 = 1
        event4s = event4s[0:nEve4]

    mjd4s = np.zeros(nEve4)
    peak4s = np.zeros(nEve4)
    rms4s = np.zeros(nEve4)
    event4s = event4s[0:nEve4]

# matches between 1 and 4
    dt14s = np.zeros(nEve1)
    ii14s = np.arange(nEve1) * 0

    if nEve4 == 1:
        nEve4 = 0

#read in all events in directory 4
    iii = 0
    for filename in event4s:
        fullname = join(dirs[3], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        mjd4s[iii] = rs.emjd
        peak4s[iii] = rs.epeak
        rms4s[iii] = rs.erms
        if (100 * int(iii/100) == iii) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))
        iii = iii + 1

    print("%5d Events in Directory: %s" % (nEve4, dirs[3]))

# there are 4 directories and 3 x 2 x 1 = 6 pairs of directories for matchs
# for each of the directories there are 4 sets of files
#   mjdNs      - full mjds of events
#   peakNs     - peaks of each events
#   rmsNs      - noise in samples near events
#   eventNs    - names of event files
# now for each of the matches there are different Indecies
# 1 = Match between 1+2
# iiMs    - index to second directory for event match
# dtMs    - time offset between first and second directory
# eventMs - name of second event that matches the first directory
# 1 = Match between 1+2
# 2 = Match between 2+3
# 3 = Match between 1+3
# 4 = Match between 1+4
# 5 = Match between 2+4
# 6 = Match between 3+4

# determine the maximum number of matches
    mjd1 = 0.
    mjd2 = 0.
    mjd3 = 0.
    mjd4 = 0.

#    maxMatch = max(nEve1, nEve2, nEve3, nEve4)+1
    maxMatch = nEve1 + nEve2 + nEve3 + nEve4 + 1
    nMatch = 0
    match1s = np.arange( maxMatch) * 0
    match2s = np.arange( maxMatch) * 0 
    match3s = np.arange( maxMatch) * 0
    match4s = np.arange( maxMatch) * 0 

    # now match event times in directories 1 and 2
    for iii in range(nEve1):
        mjd1 = mjd1s[iii]
        dtj = 0
        mjd2 = mjd2s[dtj]
        dtMin = abs(mjd1 - mjd2)
        dtPM = mjd1 - mjd2
        for jjj in range(nEve2):
            mjd2 = mjd2s[jjj]
            dt = mjd1 - mjd2
            dtabs = abs(dt)
            if dtabs < dtMin:
                dtMin = dtabs
                dtPM = dt   # dt with sign Plus/Minus = PM
                dtj = jjj
        ii12s[iii] = dtj
        dt12s[iii] = dtPM
        
    # now match event times in directories 1 and 3
    for iii in range(nEve1):
        mjd1 = mjd1s[iii]
        dtj = 0
        mjd3 = mjd3s[dtj]
        dtMin = abs(mjd1 - mjd3)
        dtPM = mjd1 - mjd3
        for jjj in range(nEve3):
            mjd3 = mjd3s[jjj]
            dt = mjd1 - mjd3
            dtabs = abs(dt)
            if dtabs < dtMin:
                dtMin = dtabs
                dtPM = dt
                dtj = jjj
        ii13s[iii] = dtj
        dt13s[iii] = dtPM

    # now match event times in directories 1 and 4
    for iii in range(nEve1):
        mjd1 = mjd1s[iii]
        dtj = 0
        mjd4 = mjd4s[dtj]
        dtMin = abs(mjd1 - mjd4)
        dtPM = mjd1 - mjd4
        for jjj in range(nEve4):
            mjd4 = mjd4s[jjj]
            if mjd4 == 0:
                dtMin = 0
                dtj = 0
                continue
            dt = mjd1 - mjd4
            dtabs = abs(dt)
            if dt < dtMin:
                dtMin = dtabs
                dtPM = dt
                dtj = jjj
        ii14s[iii] = dtj
        dt14s[iii] = dtPM

    # now match event times in directories 2 and 3
    for iii in range(nEve2):
        mjd2 = mjd2s[iii]
        dtj = 0
        mjd3 = mjd3s[dtj]
        dtMin = abs(mjd2 - mjd3)
        dtPM = mjd2 - mjd3
        for jjj in range(nEve3):
            mjd3 = mjd3s[jjj]
            dt = abs(mjd2 - mjd3)
            if dt < dtMin:
                dtMin = dt
                dtj = jjj
                dtPM = dt
        ii23s[iii] = dtj
        dt23s[iii] = dtPM

    # now match event times in directories 2 and 4
    for iii in range(nEve2):
        mjd2 = mjd2s[iii]
        dtj = 0
        mjd4 = mjd4s[0]
        dtMin = abs(mjd2 - mjd4)
        dtj = 0
        dtPM = mjd2 - mjd4
        for jjj in range(nEve4):
            mjd4 = mjd4s[jjj]
            if mjd4 == 0:
                continue
            dt = mjd2 - mjd4
            dtabs = abs(dt)
            if dtabs < dtMin:
                dtMin = dtabs
                dtPM = dt
                dtj = jjj
        ii24s[iii] = dtj
        dt24s[iii] = dtPM

    # finally match event times in directories 3 and 4
    for iii in range(nEve3):
        mjd3 = mjd3s[iii]
        dtj = 0
        mjd4 = mjd4s[dtj]
        dtMin = abs(mjd3 - mjd4)
        dtPM = mjd3 - mjd4
        for jjj in range(nEve4):
            mjd4 = mjd4s[jjj]
            dt = mjd3 - mjd4
            dtabs = abs(dt)
            if dtabs < dtMin:
                dtMin = dtabs
                dtj = jjj
                dtPM = dt
        ii34s[iii] = dtj
        dt34s[iii] = dtPM

    # now report out the minimum time offsets and times less than 1 second off
    OneMjdSec = 1./86400.
    TwoMjdSec = 2.*OneMjdSec

# now match each pair of directories.  Match 1+2
    nMatch = 0
    iMatch = 0

    for iii in range(nEve1):
        if abs(dt12s[iii]) < offset:
            dts = dt12s[iii]/OneMjdSec
            jjj = ii12s[iii]
            match1s[iMatch] = iii
            match2s[iMatch] = jjj
            iMatch = iMatch + 1
            nMatch = nMatch + 1
            if doDebug:
                print("Event %s Matches %s; Offset: %9.6f s" % (event1s[iii], event2s[jjj], dts))
                print("%5d %18.9f: %5d %18.9f" % (iii, mjd1s[iii], jjj, mjd2s[jjj]))

    nMatch12 = nMatch
    print("Found %d Event Matches between directories 1 and 2" % (nMatch12))

# find matches bewteen 1 and 3
    nMatch = 0
    for iii in range(nEve1):
        if abs(dt13s[iii]) < offset:
            dts = dt13s[iii]/OneMjdSec
            jjj = ii13s[iii]
            match1s[iMatch] = iii
            match3s[iMatch] = jjj
            iMatch = iMatch + 1
            nMatch = nMatch + 1
            if doDebug:
                print("Event %s Matches %s; Offset: %9.6f s" % (event1s[iii], event3s[jjj], dts))
                print("%5d %18.9f: %5d %18.9f" % (iii, mjd1s[iii], jjj, mjd3s[jjj]))

    nMatch13 = nMatch
    print("Found %d Event Matches between directories 1 and 3" % (nMatch13))

# find matches between 1 and 4
    nMatch = 0
    for iii in range(nEve1):
        if dt14s[iii] == 0:
            continue
        if abs(dt14s[iii]) < offset:
            dts = dt14s[iii]/OneMjdSec
            jjj = ii14s[iii]
            match1s[iMatch] = iii
            match4s[iMatch] = jjj
            iMatch = iMatch + 1
            nMatch = nMatch + 1
            if doDebug:
                print("Event %s Matches %s; Offset: %9.6f s" % (event1s[iii], event4s[jjj], dts))
                print("%5d %18.9f: %5d %18.9f" % (iii, mjd1s[iii], jjj, mjd4s[jjj]))


    nMatch14 = nMatch
    print("Found %d Event Matches between directories 1 and 4" % (nMatch14))

# find matches between 2 and 3
    nMatch = 0
    for iii in range(nEve2):
        if abs(dt23s[iii]) < offset:
            dts = dt23s[iii]/OneMjdSec
            jjj = ii23s[iii]
            match2s[iMatch] = iii
            match3s[iMatch] = jjj
            iMatch = iMatch + 1
            nMatch = nMatch + 1
            if doDebug:
                print("Event %s Matches %s; Offset: %9.6f s" % (event2s[iii], event3s[jjj], dts))
                print("%5d %18.9f: %5d %18.9f" % (iii, mjd2s[iii], jjj, mjd3s[jjj]))

    nMatch23 = nMatch
    print("Found %d Event Matches between directories 2 and 3" % (nMatch23))

# find matches between 2 and 4
    nMatch = 0
    for iii in range(nEve2):
        if dt24s[iii] == 0:
            continue
        if abs(dt24s[iii]) < offset:
            dts = dt24s[iii]/OneMjdSec
            jjj = ii24s[iii]
            match2s[iMatch] = iii
            match4s[iMatch] = jjj
            iMatch = iMatch + 1
            nMatch = nMatch + 1
            if doDebug:
                print("Event %s Matches %s; Offset: %9.6f s" % (event2s[iii], event4s[jjj], dts))
                print("%5d %18.9f: %5d %18.9f" % (iii, mjd2s[iii], jjj, mjd4s[jjj]))

    nMatch24 = nMatch
    print("Found %d Event Matches between directories 2 and 4" % (nMatch24))

# find matches between 3 and 4
    nMatch = 0
    for iii in range(nEve3):
        if abs(dt34s[iii]) < offset:
            dts = dt34s[iii]/OneMjdSec
            jjj = ii34s[iii]
            match3s[iMatch] = iii
            match4s[iMatch] = jjj
            iMatch = iMatch + 1
            nMatch = nMatch + 1
            if doDebug:
                print("Event %s Matches %s; Offset: %9.6f s" % (event3s[iii], event4s[jjj], dts))
                print("%5d %18.9f: %5d %18.9f" % (iii, mjd3s[iii], jjj, mjd4s[jjj]))

    nMatch34 = nMatch
    print("Found %d Event Matches between directories 3 and 4" % (nMatch34))

# now for all matches found, find other matches
    nMatch = iMatch
    print("Total number of pairs of matches: %d" % (nMatch))

    for iii in range(nMatch):
        mjd1 = 0
        mjd2 = 0
        mjd3 = 0
        mdj4 = 0
        if match1s[iii] != 0:
            jjj = int(match1s[iii])
            mjd1 = mjd1s[jjj] 
            event1 = event1s[jjj]
        if match2s[iii] != 0:
            jjj = int(match2s[iii])
            mjd2 = mjd2s[jjj] 
            event2 = event2s[jjj]
        if match3s[iii] != 0:
            jjj = int(match3s[iii])
            mjd3 = mjd3s[jjj] 
            event3 = event3s[jjj]
        if match4s[iii] != 0:
            jjj = int(match4s[iii])
            mjd4 = mjd4s[jjj] 
            event4 = event4s[jjj]
        avemjd = 0.5 * (mjd1+mjd2+mjd3+mjd4)

# now find new matches if not already found
        if mjd1 == 0:
            dtMin = abs(avemjd - mjd1s[0])
            dtj = 0
            for kkk in range(nEve1):
                mjd = mjd1s[kkk]
# now compute average mjd from the two found:
                dt = abs(mjd - avemjd)
                if dt < dtMin:
                    dtMin = dt
                    dtj = kkk
            if dtMin < 2.*offset:
                match1s[iii] = dtj

# now find new matches if not already found
        if mjd2 == 0:
            dtMin = abs(avemjd - mjd2s[0])
            dtj = 0
            for kkk in range(nEve2):
                mjd = mjd2s[kkk]
# now compute average mjd from the two found:
                dt = abs(mjd - avemjd)
                if dt < dtMin:
                    dtMin = dt
                    dtj = kkk
            if dtMin < 2.*offset:
                match2s[iii] = dtj

# now find new matches if not already found
        if mjd3 == 0:
            dtMin = abs(avemjd - mjd3s[0])
            dtj = 0
            for kkk in range(nEve3):
                mjd = mjd3s[kkk]
# now compute average mjd from the two found:
                dt = abs(mjd - avemjd)
                if dt < dtMin:
                    dtMin = dt
                    dtj = kkk
            if dtMin < 2.*offset:
                match3s[iii] = dtj

# now find new matches if not already found
        if mjd4 == 0:
            dtMin = abs(avemjd - mjd4s[0])
            dtj = 0
            for kkk in range(nEve4):
                mjd = mjd4s[kkk]
# now compute average mjd from the two found:
                dt = abs(mjd - avemjd)
                if dt < dtMin:
                    dtMin = dt
                    dtj = kkk
            if (dtMin < 2.*offset) and (dtj > 0):
                match4s[iii] = dtj
                print("Directory 4 has a third match: %3d %3d %3d %3d" % (match1s[iii], match2s[iii], match3s[iii], match4s[iii]))

# now finally count matches
    counts = np.arange(nMatch) * 0

    print("Finding matches")
    for iii in range(nMatch):
       count = 0
       if int(match1s[iii]) != 0:
           count = count + 1
       if int(match2s[iii]) != 0:
           count = count + 1
       if int(match3s[iii]) != 0:
           count = count + 1
       if int(match4s[iii]) != 0:
           count = count + 1
       counts[iii] = count

    countns = np.arange(4) * 0
            
    print("Finding 4 matches")
    for iii in range(nMatch):
        matchevents = [ "", "", "", "", "", "", ""]
        ematch = 0
        count = int(counts[iii])
        countns[count-1] = countns[count-1] + 1
        if count >= nPrint:
            print("Event with %d Matches" % (count))
            j1 = int(match1s[iii])
            if j1 != 0:
                print(event1s[j1], mjd1s[j1])
                matchevents[ematch] = dirs[0]+"/"+event1s[j1]
                ematch = ematch+1
            j2 = int(match2s[iii])
            if j2 != 0:
                print(event2s[j2], mjd2s[j2])
                matchevents[ematch] = dirs[1]+"/"+event2s[j2]
                ematch = ematch+1
            j3 = int(match3s[iii])
            if j3 != 0:
                print(event3s[j3], mjd3s[j3])
                matchevents[ematch] = dirs[2]+"/"+event3s[j3]
                ematch = ematch+1
            j4 = int(match4s[iii])
            if j4 != 0:
                print(event4s[j4], mjd4s[j4])
                matchevents[ematch] = dirs[3]+"/"+event4s[j4]
                ematch = ematch+1
            if doPlot:
                print("Found %d matching events: %s ... %s" % (ematch, matchevents[0], matchevents[ematch-1]))
                plotEvent = "~/bin/E %s %s %s %s" % (matchevents[0], matchevents[1], matchevents[2], matchevents[3])
                os.system(plotEvent)
                
    nSum = 0
    for iii in range(3):
        jjj = iii + 1
        print("%3d Events have %d Matches" % (countns[jjj], jjj+1))
        nSum = nSum + (jjj*countns[iii])

    print("%3d Events had no Matches" % (nEve1+nEve2+nEve3+nEve4 - nSum))

if __name__ == "__main__":
    main()

