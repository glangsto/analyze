#Python find matchs in 4 data directories
#HISTORY
#20DEC03 GIL minor updates to histogram plottoing
#20DEC02 GIL fix bugs, add histogram plotting
#20NOV30 GIL create python object to record parts of events
#20NOV29 GIL add labels, don't print duplicate matches
#20OCT05 GIL clean up code
#20OCT02 GIL count events per hour
#20AUG25 GIL look for directories up to 15
#20MAR03 GIL check matching for directories 1 and 4; code cleanup
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
import sys
import os
import time
import datetime
import glob
import numpy as np
import radioastronomy
import subprocess
import copy

# get local timezon and utc offset
batcmd="/bin/date +%Z"
timezone = subprocess.check_output(batcmd, shell=True)
parts = timezone.split()
timezone = parts[0]
#print("Zone: %s" % (timezone))

nargs = len(sys.argv)
if nargs < 2:
    print("MATCH: Match events listed in several directories")
    print("Usage: MATCH [-OF seconds] [-D] [-C Date] [dir1 dir2 dir3 dir4 ...]")
    print("Where: Optionally the user provides the maximum time offset (secs) to call a match")
    print("Where -D  Optionally print debugging info")
    print("Where -C  Optionally provide a calendar date (ie 19Nov17) instead of directories")
    print("Where -H Optionally plot histogram of events per parts of a day")
    print("Where -N <n> Optionally print matches when number is equal or greater to <n>")
    print("Where -P Plot matches")
    print("")
    print("Glen Langston, 2020 Deember 03")
    sys.exit()

# separate arguments form file names
iOffset = 1
ifile = 1
iii = ifile
offset = 1.0/86400.   # default match offset is  1 seconds = 1/86400 of a day
nday = 24             # by default divide day in 24 hours
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot
calendar = ""
doPlot = False
doHistogram = False
flagGroups = False
doDebug = False
nPrint = 4

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-OF":
        offset = np.float(sys.argv[iii+1])
        iii = iii + 1
        print("Maximum Time Offset: %8.6f s for Match: " % (offset))
        offset = offset/86400.   # convert to MJDs
        ifile = ifile + 2
    if str(anarg[0:3]) == "-N":
        nPrint = np.int(sys.argv[iii+1])
        iii = iii + 1
        if nPrint <= 0:
            nPrint = 2
        print("Print if %d or more matches" % (nPrint))
        ifile = ifile + 2
    if str(anarg[0:3]) == "-F":
        flagGroups = True
        print("Flagging groups of events")
        ifile = ifile + 1
    if str(anarg[0:3]) == "-ND":
        nday = np.int(sys.argv[iii+1])
        iii = iii + 1
        print("Divide Day into N Parts:  %d" % (nday))
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-SI":
        sigma = np.float(sys.argv[iii+1])
        iii = iii + 1
        print("Keeping Events > %7.2f Sigma " % (sigma))
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print("Note: %s" % (note))
        ifile = ifile + 2
        aFix = True
    if anarg[0:2] == "-C":
        calendar = sys.argv[iii+1]
        iii = iii + 1
        print("Matching Date: %s" % (calendar))
        ifile = ifile + 2
        aFix = True
    if anarg[0:2] == "-P":
        doPlot = True
        print("Plotting Matching events")
        ifile = ifile + 1
    if anarg[0:2] == "-H":
        doHistogram = True
        print("Plotting Histogram of events per day")
        ifile = ifile + 1
    if anarg[0:2] == "-D":
        doDebug = True
        print("Debugging")
        ifile = ifile + 1
    iii = iii + 1

nplot = 0
if nday < 1:
    nday = 24

#print("Counting Events in blocks of %5.2f hours" % (24./nday))
# 
ndirs = nargs-ifile
# 
MAXDIR = 10
NOMATCH = -9999
MAXDT = 10000.
# define the maximum number of events to keep
MAXEVENTS = 10000


def finddirs( calendar, ifile, ndirs):
    """
    finddirs() finds the names of directories matching the calendar flag
    """
    dirs = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""]
    ndir = 0
    # now find calendars
    if calendar == "":
        dirs[0] = sys.argv[ifile]
        dirs[1] = sys.argv[ifile+1]
        ndir = 2
        if n > 2:
            dir[ndir] = sys.argv[ifile+2]
            ndir = ndir+1
        if ndirs > 3:
            dir[ndir] = sys.argv[ifile+3]
            ndir = ndir+1
    else:
        for idir in range(1, 15):
            adir = "pi%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs[ndir] = adir
                ndir = ndir + 1
                print("Found directory: %s" % adir)
            if ndir >= MAXDIR:
                break

        for idir in range(1, 15):
            adir = "odroid%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs[ndir] = adir
                ndir = ndir + 1
                print("Found directory: %s" % adir)
            if ndir >= MAXDIR:
                break

        for idir in range(1, 15):
            adir = "od%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs[ndir] = adir
                ndir = ndir + 1
                print("Found directory: %s" % adir)
            if ndir >= MAXDIR:
                break
    if ndir < 2:
        print("Can not match events in less than 2 directories")
        sys.exit()
    return dirs[0:ndir]

def readEventsInDir( directory):
    """
   Read in all events in a file list and return arrays
    """

#   only examine event files
    events = list(glob.glob(os.path.join(directory,'*.eve')))
    nEve = len(events)

    if nEve < 1:
        print( "No Events in Directory; %s" % (directory))
        return nEve, "", 0., 0., 0., 0., 0.
#    else:
#        print("Found %5d events in directory: %s" % (nEve, directory))

    mjds = np.zeros(nEve)
    peaks = np.zeros(nEve)
    rmss = np.zeros(nEve)
    azs = np.zeros(nEve)
    els = np.zeros(nEve)

    rs = radioastronomy.Spectrum()
    kkk = 0
    for filename in events:
        fullname = filename
        rs.read_spec_ast(fullname)
        mjds[kkk] = rs.emjd
        peaks[kkk] = rs.epeak
        rmss[kkk] = rs.erms
        azs[kkk] = rs.telaz
        els[kkk] = rs.telel
        if (100 * int(kkk/100) == kkk) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (kkk, rs.emjd, rs.epeak, rs.erms))
        kkk = kkk + 1

    return nEve, events, mjds, peaks, rmss, azs, els
# end of readEventsInDir

# count directories with events
dirs = finddirs( calendar, ifile, ndirs)

print("Counting Events in blocks of %5.2f hours" % (24./nday))

if doDebug:
    ndir = len(dirs)
    for i in range(ndir):
        print("Dir %d: %s" % (i+1, dirs[i]))

def findpairs( EventDir1, EventDir2):
    """
    findpairs() - match events in pairs of event structures
    Inputs are Python data structures describing the events in the directories
    Outputs are indicies of matches to 2nd directory and time offsets between
    """
    # extract number of eventts in each directory
    n1 = EventDir1['n']
    n2 = EventDir2['n']
    # make sure there are plenty of indicies
    ii12s = np.zeros(n1+n2) + NOMATCH
    # init minimum time offsets
    dt12s = np.zeros(n1+n2) + MAXDT

    mjd1s = EventDir1['mjds']
    mjd2s = EventDir2['mjds']

    # now match event times in directories 1 and 2
    for lll in range(n1):
        mjd1 = mjd1s[lll]
        # assume closest match is first file in 2nd list
        dtj = 0
        mjd2 = mjd2s[dtj]
        dtMin = abs(mjd1 - mjd2)
        dtPM = mjd1 - mjd2
        for jjj in range(n2):
            mjd2 = mjd2s[jjj]
            dt = mjd1 - mjd2
            dtabs = abs(dt)
            if dtabs < dtMin:
                dtMin = dtabs
                dtPM = dt   # dt with sign Plus/Minus = PM
                dtj = jjj   # keep index to closest other event
        ii12s[lll] = dtj
        dt12s[lll] = dtPM
    return ii12s, dt12s
# end of find pairs of matches
        
def main():
    """
    Main executable for matching transient events
    """

    nDirs = len(dirs)
    nEvents = np.zeros(nDirs)
          
    iDir = 0
    # look though all directories
    for dir in dirs:
        # read all events in one directory
        nEve, events, mjds, peaks, rmss, azs, els = readEventsInDir( dir)
        # if any events
        if nEve > 0:
            # create an object with the events
            aDir = { 'dir': dir,
                     'events' : events, 
                     'n': nEve,
                     'mjds': mjds,
                     'peaks': peaks,
                     'rmss' : rmss,
                     'azs' : azs,
                     'els' : els
                     }
            if iDir == 0:
                EventDirs = { iDir : aDir }
            else:
                EventDirs.update({ iDir : aDir})
            iDir = iDir + 1
        print("%5d Events in Directory: %s" % (nEve, dir))

    # end reading all events in a directory

    # for each of the directories, count the events per hour (or other fraction of a day)
    nDir = iDir # keep track of directory count
    mjdRef = 0.
    # for each directory, find the reference mjd
    for iDir in range( nDir):
        session = copy.deepcopy(EventDirs[ iDir])
        mjds = session['mjds']
        if mjdRef == 0:
            mjdRef = int(mjds[0])
        if int(mjds[0]) < mjdRef:
            mjdRef = int(mjds[0])

    print( "MJD references: %12.6f" % (mjdRef))

    nTotal = 0
    # next, for each directory, count events in each hour/part of a day
    for iDir in range( nDir):
        session = copy.deepcopy(EventDirs[ iDir])
        nEve = session['n']
        nTotal = nTotal + nEve
        mjds = session['mjds']
        counts = np.zeros(nday)
        # now count events in each range
        for iEve in range(nEve):
            dMjd = mjds[iEve] - mjdRef
            iDay = int((dMjd*nday))
            # wrap around day fraction
            iDay = iDay % nday
            counts[iDay] = counts[iDay]+1
        session['counts'] = counts
        EventDirs[ iDir] = copy.deepcopy(session)

#       print("Hour fraction of Events in Directory %s" % (session["dir"]))

    if nDir > 0:
        session0 = EventDirs[ 0 ]
        nEve0 = EventDirs[ 0 ]['n']
        counts0 = copy.deepcopy( EventDirs[0]['counts'])
    else:
        print("No directories with events found, exiting")
        sys.exit()
    if nDir > 1:
        nEve1 = EventDirs[ 1 ]['n']
        session1 = EventDirs[ 1 ]
        counts1 = copy.deepcopy( EventDirs[1]['counts'])
    else:
        nEve1 = 0
        counts1 = np.zeros(nday)
    if nDir > 2:
        session2 = EventDirs[ 2 ]
        nEve2 = EventDirs[ 2 ]['n']
        counts2 = copy.deepcopy( EventDirs[2]['counts'])
    else:
        nEve2 = 0
        counts2 = np.zeros(nday)
    if nDir > 3:
        session3 = EventDirs[ 3 ]
        nEve3 = EventDirs[ 3 ]['n']
        counts3 = copy.deepcopy( EventDirs[3]['counts'])
    else:
        nEve3 = 0
        counts3 = np.zeros(nday)
    if nDir > 4:
        session4 = EventDirs[ 4 ]
        nEve4 = EventDirs[ 4 ]['n']
        counts4 = copy.deepcopy( EventDirs[4]['counts'])
    else:
        nEve4 = 0
        counts4 = np.zeros(nday)

    print("Hour        Telescope/Day  ")
    print("Hour      1     2     3     4")

    # prepart to compute local time 
    utcOffsetHours = time.timezone/3600. 
    utcOffsetDays = time.timezone/86400.
    utcparts = np.zeros(nday)
    print( "Time zone %s is offset %5.1f hours from UTC" % (timezone, utcOffsetHours))

    for iDay in range(nday):
        utcparts[iDay] = np.float(iDay*24./np.float(nday))
#        dayparts = utcparts - utcOffsetHours
        dayparts = utcparts
        if counts0[iDay] != 0 or counts1[iDay] != 0 or counts2[iDay] != 0 \
                or counts3[iDay] != 0 or counts4[iDay] != 0:
            if nDir == 4:
                print("%5.1f  %5d %5d %5d %5d" % (dayparts[iDay], \
                                               counts0[iDay],counts1[iDay], counts2[iDay], counts3[iDay]))
            elif nDir == 3:
                print("%5.1f %5d %5d %5d " % (dayparts[iDay], \
                                               counts0[iDay],counts1[iDay], counts2[iDay]))
            elif nDir == 2:
                print("%5.1f %5d %5d " % (dayparts[iDay], \
                                               counts0[iDay],counts1[iDay]))
            elif nDir == 1:
                print("%5.1f %5d " % (dayparts[iDay], \
                                               counts0[iDay]))
            else: # else 5 or more, just print 5
                print("%5.1f %5d %5d %5d %5d %5d" % (dayparts[iDay], counts0[iDay], \
                                                       counts1[iDay], counts2[iDay], counts3[iDay], counts4[iDay]))

# Given N telescopes, there are N!/2 pairs of observations 
# ie 2 Telescopes, 1 Pair
#    3 Telescopes, 3 Pairs
#    4 Telescopes  6 Pairs
#    5 Telescopes 30 Pairs

# determine the maximum number of matches
    maxMatch = nTotal 
    nMatch = 0

    mjd0s = EventDirs[0]['mjds']
    mjd1s = EventDirs[1]['mjds']
    mjd2s = EventDirs[2]['mjds']
    mjd3s = EventDirs[3]['mjds']

#    print( "len(mjd2s): %5d" % (len(mjd2s)))

    # for each of the pairs of directories, find closest matches
    ii01s, dt01s = findpairs( EventDirs[0], EventDirs[1])
    ii02s, dt02s = findpairs( EventDirs[0], EventDirs[2])
    ii03s, dt03s = findpairs( EventDirs[0], EventDirs[3])
    ii12s, dt12s = findpairs( EventDirs[1], EventDirs[2])
    ii13s, dt13s = findpairs( EventDirs[1], EventDirs[3])
    ii23s, dt23s = findpairs( EventDirs[2], EventDirs[3])

    # now report out the minimum time offsets and times less than 1 second off
    OneMjdSec = 1./86400.

# now match each pair of directories.  Match 1+2
# count signficant matchs
    nMatch = 0
    iMatch = 0

# offset is the minimum time, to accept as a match
# for all events in first directory, find acceptable matches in others
    for i0 in range(nEve0):
        matchcount = 1  # always match with self
        mjdave = mjd0s[i0]
        i1 = NOMATCH
        i2 = NOMATCH
        i3 = NOMATCH
        if abs(dt01s[i0]) < offset:
            i1 = int(ii01s[i0])
            mjdave = mjdave + mjd1s[i1]
            matchcount = matchcount + 1
        if abs(dt02s[i0]) < offset:
            i2 = int(ii02s[i0])
#            print("0: %5d, 1: %5d 2: %5d" % (i0, i1, i2))
            mjdave = mjdave + mjd2s[i2]
            matchcount = matchcount + 1
        if abs(dt03s[i0]) < offset:
            i3 = int(ii03s[i0])
            mjdave = mjdave + mjd3s[i3]
            matchcount = matchcount + 1
        # if any matches to this event
        if matchcount > 1:
            # a match has at least two partners
            mjdave = mjdave / np.float(matchcount)
#            if matchcount == nDir:
#                print( "MJD Ave: %12.6f; %12.6f %12.6f %12.6f, %12.6f" % (mjdave, mjd0s[i0], mjd1s[i1], mjd2s[i2], mjd3s[i3]))
            # currently can only compare 4 directories/telescopes
            match = { 'nmatch': nMatch, 'count': matchcount, 'mjd': mjdave, 'list': [ i0, i1, i2, i3] }
            if nMatch == 0:
                matchs = { nMatch: match }
            else:
                matchs.update({ nMatch : match })
            nMatch = nMatch + 1

# Now must first flag events already matched with this directory
    i0 = NOMATCH
    for i1 in range(nEve1):
        matchcount = 1
        mjdave = mjd1s[0]
        i2 = NOMATCH
        i3 = NOMATCH
        if abs(dt12s[i1]) < offset:
            i2 = int(ii12s[i1])
            mjdave = mjdave + mjd2s[i2]
            matchcount = matchcount + 1
        if abs(dt13s[i1]) < offset:
            i3 = int(ii13s[i1])
            mjdave = mjdave + mjd3s[i3]
            matchcount = matchcount + 1
        if matchcount > 1:
            match['nmatch'] = nMatch
            match['count'] = matchcount
            match['mjd'] = mjdave / np.float(matchcount)
            match['list'] = [ i0, i1, i2, i3]
            matchs[ nMatch ] = match
            nMatch = nMatch + 1

        if doDebug:
                print("Event %s Matches %s; Offset: %9.6f s" % (event1s[lll], event2s[jjj], dts))
                print("%5d %18.9f: %5d %18.9f" % (lll, mjd1s[lll], jjj, mjd2s[jjj]))
# Now must first flag events already matched with this directory
    i0 = NOMATCH
    i1 = NOMATCH
    for i2 in range(nEve2):
        matchcount = 1
        mjdave = mjd2s[0]
        i3 = NOMATCH
        if abs(dt23s[i2]) < offset:
            i3 = int(ii23s[i2])
            mjdave = mjdave + mjd3s[i3]
            matchcount = matchcount + 1
        if matchcount > 1:
            match['nmatch'] = nMatch
            match['count'] = matchcount
            match['mjd'] = mjdave / np.float(matchcount)
            match['list'] = [ i0, i1, i2, i3]
            matchs[ nMatch ] = match
            nMatch = nMatch + 1

    print("Total number of sets of matches: %d" % (nMatch))

    # now trim out all the pairs of matches that are double counted
    # ie if matched between dirs 0 1, 2 and 3
    # then must flag out matched 1,2, 1,3 and 2,3
    for lll in range((nMatch-1)):
            # now compare all pairs of matches
        matcha = matchs[ lll]
        lista = matcha['list']
        counta = matcha['count']
        # first pass, only look at 4 matches
        if counta < 4:
            continue
        for nnn in range((nMatch-lll-1)):
            kkk = lll + nnn + 1
            matchb = matchs[kkk]
            listb = matchb['list']
            countb = matchb['count']
            # count pairs of matches between two lists:
            paircount = 0
            if lista[0] == listb[0]:
                paircount = paircount+1
            if lista[1] == listb[1]:
                paircount = paircount+1
            if lista[2] == listb[2]:
                paircount = paircount+1
            if lista[3] == listb[3]:
                paircount = paircount+1
            # if the number of matchs in 2nd matches all in first
            if paircount == countb:
                countb = 0
                matchb['count'] = 0
                matchs[kkk] = matchb

    counts = np.zeros(nDir)
# finally count matches between directories
    for lll in range(nMatch):
            # now compare all pairs of matches
        matcha = matchs[ lll]
        lista = matcha['list']
        counta = matcha['count']
        if counta < 1: # if a flagged match
            continue # ignore this match
        if lista[0] != NOMATCH:
            counts[0] = counts[0] + 1
        if lista[1] != NOMATCH:
            counts[1] = counts[1] + 1
        if lista[2] != NOMATCH:
            counts[2] = counts[2] + 1
        if lista[3] != NOMATCH:
            counts[3] = counts[3] + 1

    print( " Matches of Events with other telescopes")
    print( " Tel:     2     3     4")
    for iii in range(nDir): 
        counttypes = np.zeros(nDir)
        for lll in range(nMatch):
            matcha = matchs[ lll]
            lista = matcha['list']
            counta = matcha['count'] - 1
            if counta < 0:
                continue
            if lista[iii] != NOMATCH:
                counttypes[counta] = counttypes[counta] + 1
        print( "%3d:  %5d %5d %5d" % (iii, counttypes[1], counttypes[2], counttypes[3]))

#            
    for iii in range(nDir): 
            # now compare all pairs of matches
        matcha = matchs[ lll]
        lista = matcha['list']
        counta = matcha['count']
        if counta < 1: # if a flagged match
            continue # ignore this match
        if lista[0] != NOMATCH:
            counts[0] = counts[0] + 1
        if lista[1] != NOMATCH:
            counts[1] = counts[1] + 1
        if lista[2] != NOMATCH:
            counts[2] = counts[2] + 1
        if lista[3] != NOMATCH:
            counts[3] = counts[3] + 1
#        print( "Directory: %2d - %5d"  % (iii, counts[iii]))

    nall=0
    rs = radioastronomy.Spectrum()
    # make big array to keep track of multiple matches
    files0 = EventDirs[0]['events']
    match4times = np.zeros(nMatch)
    match4gallon = np.zeros(nMatch)
    match4gallat = np.zeros(nMatch)
    match4index = np.zeros(nMatch)
    # now make a list of events that only have all matches
    for lll in range(nMatch):
        matcha = matchs[ lll ]
        lista = matcha['list']
        counta = matcha['count']
        if counta < nDir:
            continue
        # must have one directory, but others might not be present
        i0 = lista[0]
        file0 = files0[i0]
        rs.read_spec_ast(file0)
        rs.azel2radec()
        match4times[nall] = matcha['mjd']
        match4gallon[nall] = rs.gallon
        match4gallat[nall] = rs.gallat
        match4index[nall] = lll
        nall = nall + 1

    match4count = np.zeros(nall) + 1
    # compute index of coordinates within a small angular area

    dd = 0.5  # match offset is small number of degrees
    n3 = nall - 1
    # now count full, matched, events within a small angular area
    for iii in range(nall):
        # if an already counted match
        if match4count[iii] < 1:
            continue
        mmm = (nall - iii) - 1
        for jjj in range(mmm):
            kkk = iii + jjj + 1
            # if not an already flagged match
            if match4count[kkk] > 0:
#                print("Matching Matches: %d %d %d" % (iii, kkk, match4count[iii]))
                # if a match of coordinates
                if (match4gallon[kkk] - dd < match4gallon[iii] and \
                        match4gallon[kkk] + dd > match4gallon[iii] and \
                        match4gallat[kkk] - dd < match4gallat[iii] and \
                        match4gallat[kkk] + dd > match4gallat[iii]):
                    match4count[iii] = match4count[iii] + 1
                    match4count[kkk] = -1
                    match4index[kkk] = -1

    for iii in range(nall):
        # if an already counted match
        if match4count[iii] < 1:
            continue
#        print( "Event: %8.2f %5.1f %5.1f %3d" % (match4times[iii], \
#                                                      match4gallon[iii], match4gallat[iii], match4count[iii]))
    # print multiple, isolated events
    # next print complete matches within the same (small) region
    for lll in range(nMatch):
        matcha = matchs[ lll ]
        lista = matcha['list']
        counta = matcha['count']
        if counta < nDir:
            continue
        # now only print unique or 1st of many matches
        isUnique = False
        for iii in range(nall):
            if lll == match4index[iii]:
                isUnique = True
                multiples = match4count[iii]
        if not isUnique:
            continue
        # must have one directory, but others might not be present
        i0 = lista[0]
        i1 = NOMATCH
        i2 = NOMATCH
        i3 = NOMATCH
        i4 = NOMATCH
        if nDir > 1:
            i1 = lista[1]
        if nDir > 2:
            i2 = lista[2]
        i2 = lista[2]
        if nDir > 3:
            i3 = lista[3]
        if i0 == NOMATCH:
            file0 = ""
        else:
            file0 = files0[i0]
            rs.read_spec_ast(file0)
#            print("In UTC %s Ra,Dec: %7.1f,%7.1f " % (rs.utc, rs.ra, rs.dec))
#            print("In LST %8.3f %4d" % (rs.lst, i0))
            rs.azel2radec()
#            print("Ou UTC %s Ra,Dec: %7.1f,%7.1f " % (rs.utc, rs.ra, rs.dec))
#            print("Ou LST %8.3f %4d" % (rs.lst, i0))
            print("%3d 0 %s %7.1f %7.1f %7.1f %7.1f %d" %  (lll, file0, rs.ra, rs.dec, rs.gallon, rs.gallat, multiples))
            
        if i1 == NOMATCH:
            file1 = ""
        else:
            files1 = EventDirs[1]['events']
            file1 = files1[i1]
            print("%3d 1 %s" %  (lll, file1))
        if i2 == NOMATCH:
            file2 = ""
        else:
            files2 = EventDirs[2]['events']
            file2 = files2[i2]
            print("%3d 2 %s" %  (lll, file2))
        if i3 == NOMATCH:
            file3 = ""
        else:
            files3 = EventDirs[3]['events']
            file3 = files3[i3]
            print("%3d 3 %s" %  (lll, file3))
        if i4 == NOMATCH:
            file4 = ""
        else:
            files4 = EventDirs[4]['events']
            file4 = files3[i4]
            print("%3d 4 %s" %  (lll, file4))
        if doPlot:
            plotcmd = "~/Research/analyze/E %s %s %s %s" % (file0, file1, file2, file3, file4)
            os.system(plotcmd)

# now have a list of single events and multiple events
    
# now count all events happening within .1 degrees of other events.
        
# prepare to place annotations on the plot
    yoffset = np.max(counts0) + np.max(counts1) + np.max(counts2) + np.max(counts3)
# put 20 annotation lines vertically
    yoffset = yoffset / 20
    y0 = 2*yoffset
# now plot histogram of events and multiple matches
    if doHistogram:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()

# time epsilon to unhide multiple events
        adir = EventDirs[0]['dir']
        alabel = adir[0:3]
        plt.step( dayparts, counts0, where='post', label=alabel)
        max0 = np.max(counts0)
        if nDir > 1:
            adir = EventDirs[1]['dir'][0:3]
            alabel = adir[0:3]
            plt.axhline(y=max0, linestyle='dotted')
            plt.step( dayparts, counts1+max0, where='post', label=alabel)
            max0 = max0 + np.max(counts1)
        if nDir > 2:
            adir = EventDirs[2]['dir']
            alabel = adir[0:3]
            plt.axhline(y=max0, linestyle='dotted')
            plt.step( dayparts, counts2+max0, where='post', label=alabel)
            max0 = max0 + np.max(counts2)
        if nDir > 3:
            adir = EventDirs[3]['dir']
            alabel = adir[0:3]
            plt.axhline(y=max0, linestyle='dotted')
            plt.step( dayparts, counts3+max0, where='post', label=alabel)
            max0 = max0 + np.max(counts3)
        if nDir > 4:
            adir = EventDirs[4]['dir']
            alabel = adir[0:3]
            plt.axhline(y=max0, linestyle='dotted')
            plt.step( dayparts, counts4+max0, where='post', label=alabel)
            max0 = max0 + np.max(counts4)
# now annoated Az, El, Ra, Dec
        alabel = "AZ,El:  %5.1f,%5.1f" % (rs.telaz, rs.telel)
        xa = utcOffsetHours
        ya = max0 - yoffset
        ax.annotate(alabel, xy=( xa, ya), xytext=(xa , ya) )
        verticalcolors = ['r','g','b','c','m', 'k']
        # get Local Ra, Dec at UTC midnight
        utcstr = str(rs.utc)
        parts = utcstr.split()
        dateHour = "%s %5d" % (parts[0],utcOffsetHours)
        utcmidnight = datetime.datetime.strptime( dateHour, "%Y-%m-%d %H")
        print("Utc: %s, Local Midnight Utc: %s" % (utcstr, utcmidnight))

        rs.utc = utcmidnight
        rs.azel2radec()
        alabel = "RA,Dec: %5.1f,%5.1f" % (rs.ra, rs.dec)
        xa = utcOffsetHours % 24.
#        xa = utcparts[0]
        ya = max0
        ax.annotate(alabel, xy=( xa, ya), xytext=( xa, ya) )
        epsilon = 3./1440.
        # now draw vertical lines for midnight local time
        plt.axvline( utcOffsetHours, color='b', linestyle='dotted')
        alabel = "Local Midnight"
        xa = utcOffsetHours % 24.
        ax.annotate(alabel, xy=( xa, -yoffset), xytext =( xa, -yoffset))
        # now draw vertical lines for events
        for iii in range(nall):
            x4 = match4times[iii] - mjdRef
            x4 = (x4*24.) % 24.
            # now compute x position for this MJD 
            # MJD midnight = UTC midnight.   For EST, the time is actually 5 hours before midnight
#            x4 = x4 - utcOffsetHours + (epsilon*((iii%7)-3))
            x4 = x4 + (epsilon*((iii%7)-3))
            icolor = iii % 6
            plt.axvline( x4, color=verticalcolors[icolor], linestyle='dashed')
            # if a single match
            if match4count[iii] == 1:
                alabel = "%5.1f,%5.1f" % (match4gallon[iii],match4gallat[iii])
                ax.annotate(alabel, xy=( x4, y0), xytext=(x4 , y0), \
                                arrowprops=dict(facecolor='black', shrink=.003),)
            elif match4count[iii] > 1:
                alabel = "%5.1f,%5.1f (%d)" % (match4gallon[iii],match4gallat[iii], match4count[iii])
                ax.annotate(alabel, xy=( x4, y0), xytext=(x4 , y0), \
                                arrowprops=dict(facecolor='black', shrink=.003),)
            y0 = y0 + yoffset
            # keep from going off the plot
            if y0 > max0 - (2 * yoffset):
                y0 = 2 * yoffset
        plt.legend(title="Set of Obs.")

        plt.title("%s   Number of Events per hour for %d Observations" % (calendar, nDir))
        plt.xlabel("Time (UTC hours, %s offset: %2.0f hours)" % (timezone, utcOffsetHours))
        plt.ylabel("Count of Events")
        plt.show()

if __name__ == "__main__":
    main()

