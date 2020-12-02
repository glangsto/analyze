#Python find matchs in 4 data directories
#HISTORY
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
import glob
import numpy as np
import radioastronomy

nargs = len(sys.argv)
if nargs < 2:
    print("MATCH4: Match events listed in four directories")
    print("Usage: MATCH4 [-OF seconds] [-D] [-C Date] [dir1 dir2 dir3 dir4]")
    print("Where: Optionally the user provides the maximum time offset (secs) to call a match")
    print("Where -D  Optionally print debugging info")
    print("Where -C  Optionally provide a calendar date (ie 19Nov17) instead of directories")
    print("Where -N <n> Optionally print matches when number is equal or greater to <n>")
    print("Where -P Plot matches")
    print("")
    print("Glen Langston, 2020 November 30")
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
    if anarg[0:2] == "-D":
        doDebug = True
        print("Debugging")
        ifile = ifile + 1
    iii = iii + 1

nplot = 0
if nday < 1:
    nday = 24

print("Counting Events in blocks of %5.2f hours" % (24./nday))

eventCounts = np.zeros(nday)  # count events per fraction of a day
eventMatch2 = np.zeros(nday)  # count events per fraction of a day
eventMatch3 = np.zeros(nday)  # count events per fraction of a day
eventMatch4 = np.zeros(nday)  # count events per fraction of a day
eventAveAz2 = np.zeros(nday) # count events per fraction of a day
eventAveEl2 = np.zeros(nday) # count events per half hour
eventAveAz3 = np.zeros(nday) # count events per fraction of a day
eventAveEl3 = np.zeros(nday) # count events per half hour
eventAveAz4 = np.zeros(nday) # count events per fraction of a day
eventAveEl4 = np.zeros(nday) # count events per half hour
nfiles = nargs-ifile
MAXDIR = 10
NOMATCH = -9999
MAXDT = 10000.

def finddirs( calendar, ifile, nfiles):
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
        if nfiles > 3:
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
dirs = finddirs( calendar, ifile, nfiles)

print("Counting Events in blocks of %5.2f hours" % (24./nday))

if doDebug:
    ndir = len(dirs)
    for i in range(ndir):
        print("Dir %d: %s" % (i+1, dirs[i]))

def countEvents(mjd):
    """
    Update global event count for this mjd
    """
    global eventCounts

    if mjd == 0.:
        return
    iMjd = int(mjd)   # get integer day number
    iDay = mjd - iMjd # get fraction of a day
    iDay = int(iDay*nday) # get index into count array
#    print("%d %11.3f" % (iDay, mjd))
    eventCounts[iDay] = eventCounts[iDay] + 1
    return

#end of count events

def count2Matchs(mjd, az, el):
    """
    Update global event count for this mjd
    """
    global eventMatch2

    if mjd == 0.:
        return
    iMjd = int(mjd)   # get integer day number
    iDay = mjd - iMjd # get fraction of a day
    iDay = int(iDay*nday) # get index into count array
    eventMatch2[iDay] = eventMatch2[iDay] + 1
    eventAveAz2[iDay] = eventAveAz2[iDay] + az
    eventAveEl2[iDay] = eventAveEl2[iDay] + el
    return
#end of count 2 Matches

def count3Matchs(mjd, az, el):
    """
    Update global event count for this mjd
    """
    global eventMatch3

    if mjd == 0.:
        return
    iMjd = int(mjd)   # get integer day number
    iDay = mjd - iMjd # get fraction of a day
    iDay = int(iDay*nday) # get index into count array
    eventMatch3[iDay] = eventMatch3[iDay] + 1
    eventAveAz3[iDay] = eventAveAz3[iDay] + az
    eventAveEl3[iDay] = eventAveEl3[iDay] + el
    return
#end of count 3 matchs

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
        
def count4Matchs(mjd, az, el):
    """
    Update global event count for this mjd
    """
    global eventMatch4, eventAveAz4, eventAveEl4

    if mjd == 0.:
        return
    iMjd = int(mjd)   # get integer day number
    iDay = mjd - iMjd # get fraction of a day
    iDay = int(iDay*nday) # get index into count array
#    print("%d %10.2f" % (iDay, mjd))
    eventMatch4[iDay] = eventMatch4[iDay] + 1
    eventAveAz4[iDay] = eventAveAz4[iDay] + az
    eventAveEl4[iDay] = eventAveEl4[iDay] + el
    return
#end of count 4 Matchs

# define the maximum number of events to keep
MAXEVENTS = 10000

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
    nTotal = 0
    for iDir in range( nDir):
        session = EventDirs[ iDir]
        nEve = session['n']
        nTotal = nTotal + nEve
        mjds = session['mjds']
        mjd0 = mjds[0]
        counts = np.zeros(nday)
        # now count events in each range
        for iEve in range(nEve):
            dMjd = mjds[iEve] - mjd0
            iDay = int(dMjd*nday)
            # wrap around day fraction
            iDay = iDay % nday
            counts[iDay] = counts[iDay]+1
        session['counts'] = counts
        EventDirs[ iDir] = session
        print("Hour fraction of Events in Directory %s" % (session["dir"]))

    session0 = EventDirs[ 0 ]
    session1 = EventDirs[ 1 ]
    session2 = EventDirs[ 2 ]
    session3 = EventDirs[ 3 ]
    nEve0 = EventDirs[ 0 ]['n']
    nEve1 = EventDirs[ 1 ]['n']
    nEve2 = EventDirs[ 2 ]['n']
    nEve3 = EventDirs[ 3 ]['n']
    counts0 = EventDirs[0]['counts']
    counts1 = EventDirs[1]['counts']
    counts2 = EventDirs[2]['counts']
    counts3 = EventDirs[3]['counts']
    print("Hour          Telescope  ")
    print("Hour     1     2     3     4")
    for iDay in range(nday):
        if counts0[iDay] != 0 or counts1[iDay] or counts2[iDay] or counts3[iDay]:
            print("%4d %5d %5d %5d %5d" % (iDay, \
                                               counts0[iDay],counts1[iDay], counts2[iDay], counts3[iDay]))

# Given N telescopes, there are N!/2 pairs of observations 
# ie 2 Telescopes, 1 Pair
#    3 Telescopes, 3 Pairs
#    4 Telescopes  6 Pairs
#    5 Telescopes 30 Pairs

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

    maxMatch = nTotal 
    nMatch = 0
    match1s = np.arange(maxMatch) * 0
    match2s = np.arange(maxMatch) * 0 
    match3s = np.arange(maxMatch) * 0
    match4s = np.arange(maxMatch) * 0 

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
        matchcount = 0
        mjdave = mjd0s[0]
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
        if matchcount > 0:
            matchcount = matchcount + 1
            mjdave = mjdave / np.float(matchcount)
            match = { 'nmatch': nMatch, 'count': matchcount, 'mjd': mjdave, 'list': [ i0, i1, i2, i3] }
            if nMatch == 0:
                matchs = { nMatch: match }
            else:
                matchs.update({ nMatch : match })
            nMatch = nMatch + 1

# Now must first flag events already matched with this directory
    i0 = NOMATCH
    for i1 in range(nEve1):
        matchcount = 0
        mjdave = mjd1s[0]
        i2 = NOMATCH
        i3 = NOMATCH
        if abs(dt12s[i1]) < offset:
            i2 = int(ii12s[i1])
            mjdave = mjdave + mjd2s[i2]
            matchcount = matchcount + 1
        if abs(dt13s[i1]) < offset:
            i3 = int(ii13s[i1])
            mjdave = mjdave + mjd2s[i3]
            matchcount = matchcount + 1
        if matchcount > 0:
            matchcount = matchcount + 1
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
        matchcount = 0
        mjdave = mjd1s[0]
        i3 = NOMATCH
        if abs(dt23s[i2]) < offset:
            i3 = int(ii23s[i2])
            mjdave = mjdave + mjd3s[i3]
            matchcount = matchcount + 1
        if matchcount > 0:
            matchcount = matchcount + 1
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
                matchb[count] = 0
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
        print( "Directory: %2d - %5d"  % (iii, counts[iii]))

    nn=0
    for lll in range(nMatch):
        matcha = matchs[ lll]
        lista = matcha['list']
        counta = matcha['count']
        if counta < 4:
            continue
        i0 = lista[0]
        i1 = lista[1]
        i2 = lista[2]
        i3 = lista[3]

        files0 = EventDirs[0]['events']
        file0 = files0[i0]
        files1 = EventDirs[1]['events']
        file1 = files1[i1]
        files2 = EventDirs[2]['events']
        file2 = files2[i2]
        files3 = EventDirs[3]['events']
        file3 = files3[i3]
        print("%3d %s" %  (nn, file0))
        print("%3d %s" %  (nn, file1))
        print("%3d %s" %  (nn, file2))
        print("%3d %s" %  (nn, file3))
        nn = nn + 1

if __name__ == "__main__":
    main()

