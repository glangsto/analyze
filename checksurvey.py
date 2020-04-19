#Python find sets of observations making a complete map of the sky
#HISTORY
#20APR13 GIL initial version based on match 4 events
#
import matplotlib.pyplot as plt
import sys
import os
import radioastronomy
import numpy as np
import copy

nargs = len( sys.argv)
if nargs < 2:
    print("CHECKSurvey: Check for complete observations of a survey of the sky")
    print("Usage: CHECKS [-D] [-A <angle>] [-B <BeginDate>] [-E <EndDate>][dir1 dir2 ... dirN]")
    print("Where: Optionally the user provides the maximum time offset (secs) to call a match")
    print("Where -A  Optionally provide a minimum angular separation between samples (default 5 degrees)")
    print("Where -B  Optionally provide a begin calendar date (ie 19Nov17) instead of directories")
    print("Where -D  Optionally print debugging info")
    print("Where -E  Optionally provide an end calendar date (ie 19Nov22) instead of directories")
    print("Where -T 'Text' Optionally provide title for the plot")
    print("Where -P Plot Grid of gains")
    print("")
    print("Glen Langston - National Science Foundation - 2020 April 13")
    exit()


xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

# separate arguments from directory names
iOffset = 1
idir = 1
iii = idir
angle = 5.0           # min angular offset between samples on the sky
nday = 24             # by default divide day in 24 hours
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot
calendar = ""
bcalendar = ""
ecalendar = ""
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
        idir = idir + 2
    if str(anarg[0:3]) == "-N":
        nPrint = np.int( sys.argv[iii+1])
        iii = iii + 1
        if nPrint <= 0:
            nPrint = 2
        print("Print if %d or more matches" % (nPrint))
        idir = idir + 2
    if str(anarg[0:3]) == "-SI":
        sigma = np.float( sys.argv[iii+1])
        iii = iii + 1
        print("Keeping Events > %7.2f Sigma " % (sigma))
        aFix = True
        idir = idir + 2
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print("Note: ", note)
        idir = idir + 2
        aFix = True
    if anarg[0:2] == "-B":
        bcalendar = sys.argv[iii+1]
        calendar = bcalendar
        iii = iii + 1
        print("Matching Start Date: ", bcalendar)
        idir = idir + 2
        aFix = True
    if anarg[0:2] == "-E":
        ecalendar = sys.argv[iii+1]
        iii = iii + 1
        print("Matching End Date: ", ecalendar)
        idir = idir + 2
        aFix = True
    if anarg[0:2] == "-P":
        doPlot = True
        print("Plotting Matching events")
        idir = idir + 1
    if anarg[0:2] == "-D":
        doDebug = True
        print("Debugging")
        idir = idir + 1
    iii = iii + 1

nplot = 0
ndir = nargs-idir
dirs = [ "", "", "", "", "", "", "", "", "", "","", "", "", "", "", "",""]
MAXDIR = len(dirs)

nRa = int(360.01/angle)
nDec = int(180.01/angle)
print("Checking for a complete survey with samples in the RA Dec Grid: %3d,%3d" % (nRa, nDec))

# now find calendars
if bcalendar == "":
    dirs[0] = sys.argv[idir]
    if ndir > 1:
        dirs[1] = sys.argv[idir+1]
    if ndir > 2:
        dir[2] = sys.argv[idir+2]
    if ndir > 3:
        dir[3] = sys.argv[idir+3]
else:
    ndir = 0 # assume no directories, then find them
    # survey observations may be taken by a Raspberry Pi computer
    for idir in range(1,10):
        adir = "pi%d-data-" % idir
        adir = adir + bcalendar
        if os.path.isdir(adir):
            dirs[ndir] = adir
            ndir = ndir + 1
            print("Found directory: %s" % adir)
        if ndir >= MAXDIR:
            break
    # Or survey observations may be taken by an Odroid computer
    for idir in range(1,10):
        adir = "odroid%d-data-" % idir
        adir = adir + calendar
#        print("TestDir: %s" % (adir))
        if os.path.isdir(adir):
            dirs[ndir] = adir
            ndir = ndir + 1
            print("Found directory: %s" % adir)
        if ndir >= MAXDIR:
            break

    for idir in range(1,10):
        adir = "od%d-events-" % idir
        adir = adir + calendar
        if os.path.isdir(adir):
            dirs[ndir] = adir
            ndir = ndir + 1
            print("Found directory: %s" % adir)
        if ndir >= MAXDIR:
            break

if doDebug:
    for i in range(ndir):
        print("Dir %d: %s" % (i+1, dirs[i]))

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

    ndir = np.zeros(50)
    
    files1 = [f for f in listdir(dirs[0]) if isfile(join(dirs[0], f))]
    ndir[0] = len(files1)
    files2 = [f for f in listdir(dirs[1]) if isfile(join(dirs[1], f))]
    ndir[1] = len(files2)
    if len(dirs[2]) > 0:
        files3 = [f for f in listdir(dirs[2]) if isfile(join(dirs[2], f))]
        ndir[2] = len(files3)
    else:
        files3 = ""
    if len(dirs[3]) > 0: 
        files4 = [f for f in listdir(dirs[3]) if isfile(join(dirs[3], f))]
        ndir[3] = len(files4)
    else:
        files4 = ""

    count = 0
    if doDebug: 
        for i in range(4):
            print("%5d Files in Directory %d" % (ndir[i], i+1))

    nSpec1 = 0
    nSpec2 = 0
    nSpec3 = 0
    nSpec4 = 0
    nHot1 = 0
    nHot2 = 0
    nHot3 = 0
    nHot4 = 0

    rs = radioastronomy.Spectrum()

# now start reading all the files
    spec1s = files1
    spec2s = files2
    spec3s = files3
    spec4s = files4

    hot1s = files1
    hot2s = files2
    hot3s = files3
    hot4s = files4

# first directory
    for filename in files1:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "ast":
            spec1s[nSpec1] = filename
            nSpec1 = nSpec1 + 1
        if parts[nparts-1] == "hot":
            hot1s[nHot1] = filename
            nHot1 = nHot1 + 1

    if nSpec1 < 1:
        nSpec1 = 1
    mjd1s = np.zeros(nSpec1)
    peak1s = np.zeros(nSpec1)
    rms1s = np.zeros(nSpec1)
    dt12s = np.zeros(nSpec1)
    ii12s = np.arange(nSpec1) * 0
    spec1s = spec1s[0:nSpec1]
    if nSpec1 == 1:
        nSpec1 = 0

    iii = 0
    for filename in spec1s:
        fullname = join(dirs[0], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        if (100 * int(iii/100) == iii) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))
        iii = iii + 1

    if nHot1 > 0:
        hot1s = hot1s[0:nHot1]

    print("%5d Spectra in Directory: %s; Calibration Files: %5d" % (nSpec1, dirs[0], nHot1))

# second directory
    for filename in files2:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "ast":
            spec2s[nSpec2] = filename
            nSpec2 = nSpec2 + 1
        if parts[nparts-1] == "hot":
            hot2s[nSpec2] = filename
            nHot2 = nHot2 + 1

# directory 2
    if nSpec2 < 1:
        nSpec2 = 1
    mjd2s = np.zeros(nSpec2)

    if nSpec2 == 1:
        nSpec2 = 0

    if nHot2 > 0:
        hot2s = hot2s[0:nHot2]

    iii = 0
    for filename in spec2s:
        fullname = join(dirs[1], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        iii = iii + 1

    print("%5d Spectra in Directory: %s; Calibration Files: %5d" % (nSpec2, dirs[1], nHot2))

# third directory
    for filename in files3:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "ast":
            spec3s[nSpec3] = filename
            nSpec3 = nSpec3 + 1
        if parts[nparts-1] == "hot":
            hot3s[nHot3] = filename
            nHot3 = nHot3 + 1

    if nHot3 > 0:
        hot3s = hot3s[0:nHot3]

# directory 3
    if nSpec3 < 1:
        nSpec3 = 1

    for filename in spec3s:
        fullname = join(dirs[2], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        if (100 * int(iii/100) == iii) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))
        iii = iii + 1

    print("%5d Spectra in Directory: %s; Calibration Files: %5d" % (nSpec3, dirs[2], nHot3))

# forth directory
    for filename in files4:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "ast":
            spec4s[nSpec4] = filename
            nSpec4 = nSpec4 + 1
        if parts[nparts-1] == "hot":
            hot4s[nHot4] = filename
            nHot4 = nHot4 + 1

# directory 4
    if nSpec4 < 1:
        nSpec4 = 1
        spec4s = spec4s[0:nSpec4]
    if nHot4 > 0:
        hot4s = hot4s[0:nHot4]

    iii = 0
    for filename in spec4s:
        fullname = join(dirs[3], filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        if (100 * int(iii/100) == iii) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (iii, rs.emjd, rs.epeak, rs.erms))
        iii = iii + 1

    print("%5d Spectra in Directory: %s; Calibration Files: %5d" % (nSpec4, dirs[3], nHot4))

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

#    maxMatch = max(nSpec1, nSpec2, nSpec3, nSpec4)+1
    maxMatch = nSpec1 + nSpec2 + nSpec3 + nSpec4 + 1
    nMatch = 0
    match1s = np.arange( maxMatch) * 0
    match2s = np.arange( maxMatch) * 0 
    match3s = np.arange( maxMatch) * 0
    match4s = np.arange( maxMatch) * 0 


if __name__ == "__main__":
    main()

