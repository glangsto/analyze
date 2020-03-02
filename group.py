#Python Script to sort events into long and short duration events
#HISTORY
#19NOV07 GIL determine the durations of individual groups
#19OCT21 GIL inital version based on matchevents.py
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
offset = 240.         # default group match is 4 minutes
ngroup = 4            # by default 4 or more events is a group
nblock = 256          # number of samples to FFT
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot
doDebug = False
doPlot = False
doMove = False

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-OF":
        offset = np.float( sys.argv[iii+1])
        iii = iii + 1
        print("Time Offset Interval: %8.6f s for Group: " % ( offset))
        ifile = ifile + 2
    if str(anarg[0:3]) == "-G":
        ngroup = np.int( sys.argv[iii+1])
        iii = iii + 1
        if ngroup < 1:
            print("Minimum group size is 2")
            ngroup = 2
        print("N of events to consider a long group: ", ngroup)
        ifile = ifile + 2
    if str(anarg[0:3]) == "-SI":
        sigma = np.float( sys.argv[iii+1])
        iii = iii + 1
        print("Keeping Events > %7.2f Sigma " % (sigma))
        ifile = ifile + 2
    if str(anarg[0:2]) == "-M":
        print("Moving Events")
        doMove = True
        ifile = ifile + 1
    if str(anarg[0:3]) == "-D":
        print("Debugging")
        doDebug = True
        ifile = ifile + 1
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print("Note: ", note)
        ifile = ifile + 2
    if anarg[0:2] == "-P":
        doPlot = True
        print("Plotting Ungroupded Events")
        ifile = ifile + 1
    iii = iii + 1

N = nblock                 # abreviation
nplot = 0

MAXEVENTS = 5000

def moveFiles( inDir, groupDir, nInGroup, groupFiles):
    """
    Move files to the group subdirectory
    """
    global doDebug
    outDir = "./" + groupDir
    if doDebug:
        print("Moving %d events to directory %s " % ( nInGroup, groupDir))
    mkdir = "mkdir %s 2> /dev/null" % (outDir)
    if not os.path.isfile( outDir):
        os.system(mkdir)

    for iii in range(nInGroup):
        inFile = "./" + inDir + "/" + groupFiles[iii]
        mvEvent = "mv %s %s 2> /dev/null" % (inFile, outDir)
        if doDebug:
            # print mvEvent
            doDebug = True
        else:
            os.system(mvEvent)

    return

def main():
    """
    Main executable for gridding astronomical data
    """
    global offset

    nargs = len(sys.argv)
    if nargs < 2:
        print('GROUP: GROUP events into long and short duration')
        print('usage: GROUP [-OF seconds] [-G number] [-M] [-P] [-D] dir1')
        print('where: -M  move groups into a sub-directory')
        print('       -P  plot ungrouped events')
        print('       -G  <number> is the minimum number of events to call a group')
        print('       -D  print debugging info')
        exit()

    dir1 = sys.argv[ifile]
    print("Dir 1: ", dir1)

    from os import listdir
    from os.path import isfile, join
    files1 = [f for f in listdir(dir1) if isfile(join(dir1, f))]
    
    groupDir = dir1 + "/groups"

    print("Time offset to identify a group: %7.2f s " % (offset))

    offset = offset/86400.   # convert to MJDs

    print("%5d Files in Directory 1" % (len(files1)))
#    print files
    nGroups = 0

    rs = radioastronomy.Spectrum()

    event1s = files1
    nEve1 = 0

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

    
    # now extract only the list of events, sorted in time order
    event1s = sorted( event1s[0:nEve1])
    mjd1s = np.zeros(nEve1)
    peak1s = np.zeros(nEve1)
    rms1s = np.zeros(nEve1)

    # now get event signficances
    iii = 0
    for filename in event1s:
        fullname = join(dir1, filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        # if this is a significant event
        if rs.epeak / rs.erms > sigma:
            mjd1s[iii] = rs.emjd
            peak1s[iii] = rs.epeak
            rms1s[iii] = rs.erms
            event1s[iii] = filename 
            iii = iii + 1
        if 50 * int(iii/50) == iii:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))

    print("Read %5d events, of which %5d are > %7.2f sigma" % (nEve1, iii, sigma))

    # new count of significant events
    nEve1 = iii

    inGroup = False
    nInGroup = 0
    nAlone = 0
    groupFiles = copy.deepcopy( event1s)
    aloneEvents = copy.deepcopy( event1s)


    # now group event times

    lastt = mjd1s[0]

    if doDebug:
        print(" Offset limit for grouping: %7.2fs" % ( offset*86400.))

    iFirst = 0
    iLast = 0

    for iii in range(nEve1):
        t = mjd1s[iii]
        dt = abs(t - lastt)

        if doDebug:
            print("Event %d offset from previous: %7.2f (%d): " % ( iii, dt * 86400.,  nInGroup))
        if dt == 0.:
            inGroup = False
            groupFiles[0] = event1s[iii]
            nInGroup = 1
            firstFile = event1s[iii]
            iFirst = iii
            lastFile = event1s[iii]
            iLast = iii
        elif dt < offset:
            inGroup = True
            groupFiles[nInGroup] = event1s[iii]
            nInGroup = nInGroup + 1
            lastFile = event1s[iii]
            iLast = iii
        else:
            if doDebug:
                print("N in Group: %d" % (nInGroup))
            # if was in a group, but no longer
            if nInGroup >= ngroup:
                if doMove:
                    moveFiles( dir1, groupDir, nInGroup, groupFiles)
                dMjd = mjd1s[iLast] - mjd1s[iFirst]
                duration = 86400. * dMjd
                print("Group of %5d events (%s-%s), dur: %9.3f (s)" % (nInGroup, firstFile, lastFile, duration))
                nGroups = nGroups + 1
            else:
                kkk = iFirst
                if doDebug:
                    print(" First: %d, Last: %d"  % (iFirst, iLast))
                for jjj in range(nInGroup):
                    if doDebug:
                        print(" Index: %d, nAlone: %d; File Index %d" % (jjj, nAlone, kkk))
                    aloneEvents[nAlone] = event1s[kkk]
                    nAlone = nAlone + 1
                    kkk = kkk + 1
            # start a new group with this file
            groupFiles[0] = event1s[iii]
            nInGroup = 1
            firstFile = event1s[iii]
            iFirst = iii
            lastFile = event1s[iii]
            iLast = iii
            inGroup = False
        lastt = t
            
    # if was in a group, but no longer
    if nInGroup >= ngroup:
        if doMove:
            moveFiles( dir1, groupDir, nInGroup, groupFiles)
        nGroups = nGroups + 1
    else:
        kkk = iFirst
        for jjj in range(nInGroup):
            aloneEvents[nAlone] = event1s[kkk]
            nAlone = nAlone + 1
            kkk = kkk + 1

    print("Found %5d groups and %5d isolated events " % (nGroups, nAlone))
    return

if __name__ == "__main__":
    main()

