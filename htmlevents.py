#Python find matchs in data directories and write html
#HISTORY
#25Dec15 GIL separate out plotHistograms and update for new arguments
#25Dec12 GIL copy matched event files to directories
#25Dec13 GIL write histograms to directory of date
#25Dec12 GIL copy matched event files to directories
#25Dec11 GIL create HTML output logs.
#25Dec02 GIL Sun RA and Dec postion match telescope
#25Nov19 GIL fix mjds index out of range
#25Oct06 GIL log successful fits to event groups, use getUtcOffset
#25Oct03 GIL fit group matches
#25Sep30 GIL group matches
#25Sep29 GIL group events and catagorize
#25Sep24 GIL modularize
#25Sep23 GIL find indexing problem with matched events.
#25Sep17 GIL check match lists for repeats
#25Sep16 GIL finish cleanup and repeat matching code
#25Sep12 GIL Do not match one event to many events in other telescopes.
#25Sep05 GIL optionally fit gaussians to event data streams
#25Aug21 GIL Fix matching, input to findTime was backwards
#25Aug20 GIL Modularize a bit, move findDate to separate file
#25Aug11 GIL test for single events, to avoid len(1) failing
#25Aug11 GIL test for single events, to avoid len(1) failing
#25Aug10 GIL speed up matching of events, using better sorter
#25Aug09 GIL use file name for event time, instead of reading files.
#25Jul06 GIL write outputs to OUTPUT directory
#25May02 GIL add more x ticks and sub-ticks
#25Mar20 GIL add Verbose (-V) to print each file being examined
#25Feb22 GIL check for events in directories 10-15
#24Oct03 GIL add more galactic annotations
#24Sep03 GIL log matchs
#24May06 GIL enable limiti9ng matches due to sigma
#24Jan26 GIL print output file pdf
#23Jul26 GIL compute average RA, Dec for groups of events
#23Jul25 GIL add 1 hour x tick interval
#23Jun07 GIL enable flagging transit of a specific location
#23Mar30 GIL change histogram plot order to have label match plot order
#23Mar21 GIL fix plotting of more than 4 directories
#23Mar03 GIL Fix matching with only 2 directories
#23Feb24 GIL minor plotting improvements
#23Feb15 GIL numpy revisions
#22May20 GIL fix MATCH plotting the 5th telecope, fix UTC offset to timezone
#22Apr07 GIL fix histogram alignment
#21Dec14 GIL color in the bars in the counts/hour plot
#21Sep01 GIL ignore data for El < 0
#21MAY10 GIL clean up sun plotting
#21APR30 GIL fix sun azimuth when pointed north
#21APR28 GIL fix case of absolutely no matches
#21APR21 GIL force azimuth = 180 for Sun position calculations
#21MAR21 GIL Increase interval between Galactic coordinate calc.
#21MAR19 GIL Complete code for 5 directories/telescopes
#21FEB26 GIL Show Galactic Plane range
#21FEB25 GIL FInd the SUN transit time
#21FEB24 GIL show LST of SUN
#21FEB23 GIL remove extra character in time zone
#20DEC04 GIL allow histogram plotting in case of only one directory
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
import datetime
import copy
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.coordinates import AltAz
import radioastronomy
from jdutil import date_to_mjd
from findDate import findDirs, readEventsInDir
from findMatches import *
from groupMatches import *
import utcOffset
from plotHistogram import *
from eventDay import eventDaySummary

nargs = len(sys.argv)
if nargs < 2:
    print("HTMLMATCH: Match events listed in several directories")
    print("Usage: HTMLMATCH [-OF seconds] [-D] [-C Date] [-E] [-N <n>] [-G] [dir1 dir2 dir3 dir4 ...]")
    print("Where:")
    print(" -OF Optionally the user provides the maximum time offset (secs) to call a match")
    print(" -C  Optionally provide a calendar date (ie 19Nov17) instead of directories")
    print(" -D  Optionally print debugging info")
#    print(" -E Optionally only show Events with elevation above zero")
    print(" -F Optionally Do Not flag groups of events")
    print(" -G Optionally plot +/- 10 degrees of galactic plane")
    print(" -H Optionally plot histogram of events per parts of a day")
    print(" -L <logFileName> Set directory/file for event matches")
    print(" -M <ra> <dec> Optionally mark the transit of a coordiante")
    print(" -N <n> Optionally print matches when number is equal or greater to <n>")
    print(" -O <directory> write plots to output directory")
    print(" -ND <n> Optionally divide the day into <n> parts, default 24")
    print(" -P Plot matches")
#    print(" -SI <min sigma> Only match if event obove <min sigma>")
    print(" -Q Quiet processing, no showing plots, only save .pdf, .svg")
    print(" -Y <offset> Optionally adds an offset to the event plots (-P option)")
    print("")
    print("")
    print("Matching of events from several telescopes is more complex")
    print("than you might think.  It is fairly easy to find the closest")
    print("Events in other telescopes to a single event in one telescope.")
    print("But it is assumed")
    print("that if there are many matching events in other telescopes, only")
    print("the closest in time for each telescope is the match, not all close")
    print("Events.  Filtering requires finding the closest event for each")
    print("Telescope, then counting and removing the other close matches for")
    print("the second telescope")
    print("")
    print("Glen Langston, 2025 September 12")
    sys.exit()

# separate arguments from file names
ifile = 1
iii = ifile
OneMjdSec = float( 1./86400.)
tOffset = OneMjdSec   # default match offset is 1/2 seconds = 1/86400 of a day
nDay = 24             # by default divide day in 24 hours
sigma = 4.0           # Minimum sigma to accept
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot
calendar = ""
doPlot = False
doPlot = True
doQuiet = False
doHistogram = False
doLog = True
logFileName = ""
outDirName = "~/match"
flagGroups = False
yoffset = 0.0
doGalactic = True
doGaussian = False
#optinally mark transit of a known object of interest
doTransit = False
transitObject = ""
raTransit = 0.0     # ra and dec in degrees of object
decTransit = 0.0
# Limit diagnostic print output
nPrint = 4
minEl = -100.
# temporarily turn on/off debugging
verbose = True
verbose = False

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-O":  # if output file provided
        iii = iii + 1
        outDirName = str(sys.argv[iii])
        print( "Writing output to %s" % (outDirName))
    if str(anarg[0:3]) == "-OF":
        tOffset = float(sys.argv[iii+1])
        iii = iii + 1
        print("Maximum Time Offset: %8.6f s for Match: " % (tOffset))
        tOffset = tOffset/86400.   # convert to MJDs
        ifile = ifile + 2
    if str(anarg[0:3]) == "-N":
        nPrint = int(sys.argv[iii+1])
        iii = iii + 1
        if nPrint <= 0:
            nPrint = 2
        print("Print if %d or more matches" % (nPrint))
        ifile = ifile + 2
    if str(anarg[0:3]) == "-E":
        minEl = 1.
        print("Only counting events when telescope > %5.1f (d) el." % (minEl))
        ifile = ifile + 1
    if str(anarg[0:3]) == "-F":
        flagGroups = True
        print("Not Flagging groups of events")
        ifile = ifile + 1
    if str(anarg[0:3]) == "-G":
        doGaussian = True
        print("Fitting Gaussians to event sequences")
        ifile = ifile + 1
    if str(anarg[0:3]) == "-L":  # if logging events
        iii = iii + 1
        logFileName = sys.argv[iii]
        print("Logging Events in: %s" % (logFileName))
        doLog = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-ND":
        iii = iii + 1
        nDay = int(sys.argv[iii])
        print("Divide Day into N Parts:  %d" % (nDay))
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:2]) == "-M":   # if marking transit of a location
        doTransit = True
        iii = iii + 1
        transitObject =  sys.argv[iii]  # name of object of interest
        iii = iii + 1
        raTransit = float(sys.argv[iii])
        iii = iii + 1
        decTransit = float(sys.argv[iii])
        print("Marking Transit of %s, RA, Dec: %7.2f,%7.2f (deg)" % \
              (transitObject, raTransit, decTransit))
        ifile = ifile + 3
    if str(anarg[0:3]) == "-SI":
        sigma = float(sys.argv[iii+1])
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
    if anarg[0:2] == "-Q":
        doQuiet = True
        print("Quiet running, without showing plots")
        ifile = ifile + 1
    if anarg[0:2] == "-H":
        doHistogram = True
        print("Plotting Histogram of events per day")
        ifile = ifile + 1
    if anarg[0:2] == "-D":
        doDebug = True
        print("Debugging")
        ifile = ifile + 1
    if anarg[0:2] == "-V":
        verbose = True
        print("Verbose output")
        ifile = ifile + 1
    if anarg[0:2] == "-Y":
        iii = iii + 1
        yoffset = float(sys.argv[iii])   #plot offset for events
        ifile = ifile + 2
    iii = iii + 1

eventLogName = "~/match/EventMatches.log"
nplot = 0
# if dividing day into less that 1 part, divide day into hours
if nDay < 1:
    nDay = 24
# number of directories provided
nDir = nargs-ifile
#
MAXDIR = 20
NOMATCH = -9999
MAXDT = 365.
# Only read one event file, to get location and az, el
readAzEl = False

 
def mjd_to_ymd(mjd):
    """
    Converts a Modified Julian Date (MJD) to a year, month, and day tuple.
    """
    # MJD 0 corresponds to November 17, 1858, 00:00:00 UTC
    mjd_epoch = datetime.datetime(1858, 11, 17, 0, 0, 0)

    # Calculate the timedelta from the MJD epoch
    # We subtract 0.5 from MJD because MJD starts at midnight, while JD starts at noon.
    # However, since we are directly working with datetime objects and their difference in days,
    # we can simply add the MJD value as days to the epoch.
    target_date = mjd_epoch + datetime.timedelta(days=mjd)

    return target_date.year, target_date.month, target_date.day, \
        target_date.hour, target_date.minute, target_date.second

def logEvent( logFile, files, event, gallon, gallat):
    """
    logEvent() adds one event to the event log
    where
    logFile = file pointer to open log file
    files   = composed of this event
    event   = event to be added to the log
    """

    iEvent = int( event['nmatch'])
    aFlag  = event['flag']
    nTel   = int( event['count'])
    print("#Event: %5d nTel: %2d flag: '%s' Gal Lon,lat: %7.2f %7.2f" % ( iEvent, nTel, aFlag, gallon, gallat), file=logFile)
    for afile in files:
        # aline file names by adding a space for shorter names
        fileparts = afile.split("-")
        aTel = fileparts[0]
        if len(aTel) < 4:
            aspace = " "
        else:
            aspace = ""
        print("%s%s" % (aspace, afile), file=logFile)

## return end of logEvent()

def logGroups( logFile, calendar, nDir, nGroup, nTen, nHundred, nThousand):
    """
    logGroups() = writes summary of the events and groups found for this date
    where
    logFile = pointer to open logr file
    nGroup  = total number of groups found
    nTen    = number of groups with between 10 and 99 members   (X)
    nHundred = number of groups with 100 to 999 members         (C)
    nThousand = number of groups with more than 999 members     (M)
    """
    print( "", file=logFile)
    print(     "#DATE    =%s" % (calendar), file=logFile)
    print(     "#NTEL    =%5d" % (nDir), file=logFile)
    print(     "#NMATCH  =%5d" % (nGroup), file=logFile)
    if nTen > 0:
        print( "#TEN     =%5d" % (nTen), file=logFile)
    if nHundred > 0:
        print( "#HUNDRED =%5d" % (nHundred), file=logFile)
    if nThousand > 0:
        print( "#THOUSAND=%5d" % (nThousand), file=logFile)
#    print( "#  Tel   Directory             Event File            Ra     Dec     GLon    GLat N Flashes", file=logFile)

    return
    # end of logGroups()

def htmlGroups( directoryDate, calendar, nDir, \
                nGroup, nTen, nHundred, nThousand, \
                nFit, fits, mjdRef):
    """
    htmlGroups() creates a web page for this days events
    where
    directoryDate = directory to contain web page
    nGroup  = total number of groups found
    nTen    = number of groups with between 10 and 99 members   (X)
    nHundred = number of groups with 100 to 999 members         (C)
    nThousand = number of groups with more than 999 members     (M)
    n
    """

    htmlName = "%s/%s.html" % (directoryDate, calendar)
    htmlFile = open( htmlName, 'w')

    print("<html>Summary of events on %s </html>" % (calendar), file=htmlFile)
    print("<h1>Summary of events detected by %d telescopes on %s<h1>" %
          (nDir, calendar), file=htmlFile)

    if nGroup <= 0:
        print("<p><h2>No Events matched on %s</h2>" % (calendar),
          file=htmlFile)
    else:
        print("<p><h2>%d Events matched on %s</h2>" % (nGroup, calendar),
          file=htmlFile)

        print("<p><table border='2'>",
              file=htmlFile)
        print("<tr><td>Day</td><td>Isolated</td><td>Few</td><td>Groups</td><td>Groups</td><td>Groups of</td></tr>",
              file=htmlFile)
        print("<tr><td>    </td><td>Events</td>  <td>Events</td><td>10 or more</td><td>100 or more</td><td> 1000 or more</td></tr>",
              file=htmlFile)

        print("<tr><td> %s</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td></tr>" % (calendar, nGroup, nTen, nHundred, nThousand),
              file=logFile)
        # finish table
        print("</table><p>", file=htmlFile)

    print('<a href="../index.html"> Year </a>', file=htmlFIle)
    print( "", file=logFile)
    close( htmlFile)
    
    return
    # end of htmlGroups()

def logFits( logFile, nFit, fits, mjdRef):
    """
    logFits() adds successful fits to groups of events to the event log
    where
    logFile  = open file pointer to log file
    rs       = radio spectrum of one of the files for the fit
    nFit     = number of succesful fits
    fits     = array of fit objects
    mjdRef   = reference MJD for this observing session
    """
    print( "#NFIT    =%5d" % (nFit), file=logFile)
    for iFit in range(nFit):
        time = fits[iFit]['time']
        mjd = mjdRef + time

        # Convert MJD to astropy Time object (assumes UTC scale)
        time_obj = Time(mjd, format='mjd', scale='utc')

        # Convert astropy Time object to Python datetime object (in UTC)
        utc_datetime = time_obj.to_datetime()

        # Convert datetime object to a floating-point Unix timestamp (seconds since epoch)
        utc_float_timestamp = utc_datetime.timestamp()
        utcYYmmmDD = utcOffset.strAst( utc_datetime)

        utcstr = str(utc_datetime)
        print(f"MJD: {mjd}")
        print(f"UTC Datetime: {utc_datetime}")
        print(f"UTC AST time: {utcYYmmmDD}")
        print(f"UTC Floating-Point Timestamp: {utc_float_timestamp}")

        peak  = fits[iFit]['peak']
        rms   = fits[iFit]['rms']
        stdDev  = fits[iFit]['stdDev']
        fwhm  = fits[iFit]['fwhm']
        print("Fit%2d %7.1f+/-%4.0f %s" % (iFit, peak, rms, utcYYmmmDD), file=logFile)
        print("Time  %7.1f+/-%5.1f hours == +/- %6.2f minutes" % (time*24., fwhm*24., fwhm*1440.), file=logFile)

    return
    # end of logFits()

# count directories with events
inList = sys.argv[ifile:-1]
# count number of directories events matching input

if nDay != 24:
    nminutes = int(1440.001/nDay)
    print("Counting Events in blocks of %5d minutes" % (nminutes))

def matchClose( nall, matchcount, matchindex, dd = 10.*60./86400.) :
    """
    Count nearby in time, but not close enough to be an exact match
    """
    # for all events,
    for iii in range(nall):
        # if an already counted match
        if matchcount[iii] < 1:
            continue
        mmm = (nall - iii) - 1
        # look at all other events
        for jjj in range(mmm):
            # starting at one after the current event
            kkk = iii + jjj + 1
            # if not already counted as an event
            if matchcount[kkk] > 0:
                # if a match of time
                matchcount[kkk] = -1
                matchindex[kkk] = -1
    return matchcount, matchindex
    # end of matchClose()

def writeHornSummary( summaryFile, dirname, nShort, nLong, az, el):
    """
    write a summary of the events detected by one horn telescope
    inputs:
    summaryfile - an open file for write with append
    dirname - telescope event directory, includes name and date
    nShort - number of short events seen
    nLong - number of long events seen
    """
    outstr = "%s %5d %5d %6.1f %6.1f" % (dirname, nShort, nLong, az, el)

    summaryFile.write(outstr)
    # end of write file
    return

def main():
    """
    Main executable for matching transient events
    """

# dirs is the list of names of directories that match the input string.
# nDir is the count of directory names
    dirs, nDir = findDirs( inList, calendar)
    nEvents = np.zeros(nDir, dtype=int) # init the counts of events in each directory

    iDir = 0
    # look though all directories
    for dir in dirs:
        # read count of events, event file names and MJDs of events
        nEve, eventNames, mjds = readEventsInDir( dir)
        # convert List to array
        mjds = np.array( mjds)
        # if any events
        if nEve > 0:
            # create an object with the events
            aDir = { 'dir': dir,              # name of directory
                     'events' : eventNames,   # names of events in directories
                     'n': nEve,               # count of events in directory
                     'mjds': mjds,            # array of MJDs of these events
                     }
            fullName = dir[0:4]
            fullParts = fullName.split("-")
            telName = fullParts[0]
            if iDir == 0:                     # initialize a list of directories
                eventDirs = { iDir : aDir }
                telNames = [ telName ]
            else:
                eventDirs.update({ iDir : aDir})   #else add to directory list
                telNames.append( telName)
            iDir = iDir + 1
            print("%5d Events in Directory: %s (%s)" % (nEve, dir, telName))

    # end reading all events in a directory

    # for each of the directories, count the events per fraction of a day)
    nDir = iDir # keep track of directory count
    mjdRef = 0.
    # for each directory, find the reference mjd, which is earliest mjd in list
    for iDir in range( nDir):
        mjds =  eventDirs[ iDir]['mjds']
        if isinstance( mjds, np.ndarray):
            nmjds = len(mjds)
            anmjd = 0.
            if iDir < nmjds:
                anmjd = int(mjds[iDir])
        else:    
            anmjd = mjds
        if mjdRef == 0:        # if no MJDs yet, use first
            mjdRef = anmjd
        if anmjd < mjdRef:     # if this MJD is earlier, use it
            mjdRef = anmjd

    if nDir > 0:
        print( "MJD reference: %12.6f" % (mjdRef))

    rs = radioastronomy.Spectrum()
    nTotal = 0   # count total number of events in all directories
    # next, for each directory, count eventNames in each hour/part of a day
    totalCountsInDt = np.zeros(nDay, dtype=int)   # initialize the counts per fraction of a day
    for iDir in range( nDir):
        filename = eventDirs[ iDir]['events'][0]
        rs.read_spec_ast( filename)
        # examine each directories events.  A session is one telescope-day
        session = copy.deepcopy(eventDirs[ iDir])
        nEve = session['n']
        nTotal = nTotal + nEve
        mjds = session['mjds']
        countsInDt = np.zeros(nDay, dtype=int)   # initialize the counts per fraction of a day
        # now count events in each range
        for iEve in range(nEve):
            dMjd = mjds[iEve] - mjdRef
            # compute index into counts as a fraction of a day.
            iDay = int((dMjd*nDay))
            # wrap around day fraction
            iDay = iDay % nDay
            # these counts are only for a single telescope
            # to plot events as a function of time.
            countsInDt[iDay] = countsInDt[iDay]+1
        # a session is one day and one telescope
        session['counts'] = countsInDt
        eventDirs[ iDir] = copy.deepcopy(session)

    if nDir > 0:                 # deal with finding no telescopes
        print("Count of Events in Time interval for each Telescope ")
        print("                    Telescope")
        print(" Time ", end="")

    nEves = np.zeros(nDir)
    # now for all time intervals show events in that interval
    for iDir in range( nDir):
        nEves[iDir] = eventDirs[ iDir ]['n']
        countsInDt = copy.deepcopy( eventDirs[iDir]['counts'])

        if iDir == 0:
            countsMatrix = [ countsInDt ]
        else:
            countsMatrix.append( countsInDt)
        print("%5s " % (telNames[iDir]), end = "")
    if nDir > 0:
        print("")

    # init total per telescope
    telTotals = np.zeros(nDir)
    # now print events as a function of time, for each telescope if any has an event
    for iDay in range(nDay):
        eventsFound = False
        for iDir in range(nDir):
            if countsMatrix[iDir][iDay] != 0:
                eventsFound = True
                continue
        if eventsFound:
            utcN = float(iDay*24./float(nDay))
            print("%5.1f " % (utcN), end = "")
            for iDir in range(nDir):
                print("%5d " % (countsMatrix[iDir][iDay]), end =  "")
                telTotals[iDir] = telTotals[iDir] + countsMatrix[iDir][iDay]
                totalCountsInDt[iDay] = totalCountsInDt[iDay] + countsMatrix[iDir][iDay]
            print("")

    print("Total:", end = "")
    for iDir in range( nDir):
        print("%5d " % (int(telTotals[iDir])), end = "")
    if iDir < 1:
        print("    0")
    else:
        print("")
        # there are (nDir-1)! possible pairs
        print("       Count of matches with other telescopes")

    if nDir < 1:
        exit()
    nMatch, eventList = findMatches( nDir, eventDirs, telNames, tOffset = tOffset, verbose=False)

    nMin = 3
    nRemain, eventFinal = compressEvents( nMatch, eventList, nMin=nMin, verbose = False)
    print("##############################################################################")
    print("Compressed %d events to %d events with minimum number of telescopes=%d" % \
          (nMatch, nRemain, nMin))

    print("##############################################################################")
    eventTrim = trimEvents( nRemain, eventFinal, nDir, eventDirs)

    verbose = False

    if verbose:
        for iGroup in range( nRemain):
            aEvent = eventTrim[iGroup]
            aveMjd = aEvent['mjd']
            matches = aEvent['list']
            showMatch( iGroup, aveMjd, nDir, matches, eventDirs)

    # group matches before plotting
    print("##############################################################################")
    nGroup, nTen, nHundred, nThousand, groupEvents = groupMatches( \
                        mjdRef, nRemain, nDir, eventTrim, nDay, verbose = True)
    print("Grouped %d events into %d Groups" % (nRemain, nGroup))
    print("Groups, %3d with > 10 members, %d with > 100 and %d > 1000 members" % \
          (nTen, nHundred, nThousand))

    # create full path and open log file
    logFile, directoryDate = logDirectory( calendar, outDirName, \
                                           eventLogName)
    if doLog:
        logGroups( logFile, calendar, nDir, nGroup, nTen, nHundred, nThousand)

    rs.azel2radec()

    if nTen + nHundred + nThousand > 0:
        nFit, fits = groupFit( nDay, totalCountsInDt, nGroup, groupEvents, nFitMax = 3, verbose=True)
        if nFit > 0:
            print(" Found %3d fits to groups of events" % (nFit))
            logFits( logFile, nFit, fits, mjdRef)

    matchtimes = np.zeros(nGroup)
    matchcounts = np.zeros(nGroup)
    matchgallon = np.zeros(nGroup)
    matchgallat = np.zeros(nGroup)

    # for the time being force debug
    verbose = True

    # initialize the file list
    files = []
    groupFlags = [ ]
    # now fill arrays with coordinates and info, for plotting
    for iGroup in range(nGroup):
        aGroup = groupEvents[ iGroup]
        groupFlags.append( aGroup['flag'])
        matches = aGroup['list']
        # count becomes number of events in a group
        nEvents = int(aGroup['count'])
        matchtimes[iGroup] = aGroup['mjd']
        matchcounts[iGroup] = nEvents
        files = []
        fileNames = " " 
        nTel = 0
        aveGalLon = 0.
        aveGalLat = 0.
        for iDir in range(nDir):
            iMatch = matches[iDir]
            if iMatch != NOMATCH:
                # find file name in sets of names
                fileName = eventDirs[iDir]['events'][iMatch]
                if doPlot or doLog:
                    files.append( fileName)
                    fileNames = fileNames + " " + fileName
                rs.read_spec_ast( fileName)
                rs.azel2radec()
                aveGalLon = aveGalLon + rs.gallon
                aveGalLat = aveGalLat + rs.gallat
                nTel = nTel + 1
        # now pass info
        if nTel > 0:
            matchgallon[iGroup] = aveGalLon/nTel
            matchgallat[iGroup] = aveGalLat/nTel
        else:
            print("No Telescopes for Event/Group %d!" % (iGroup))
        aGroup['count'] = nTel

        if verbose:
            aveMjd = aGroup['mjd']
            showMatch( iGroup, aveMjd, nDir, matches, eventDirs)

        # now copy matched file to shared directory
        copyMatch( iGroup, directoryDate, aveMjd, nDir, matches, eventDirs)

        if doLog:
              logEvent( logFile, files, aGroup, \
                        matchgallon[iGroup], matchgallat[iGroup])
        if doPlot:
            plotcmd = "~/Research/analyze/E -O %s -Y %.2f %s " % \
                (directoryDate, yoffset, fileNames)
            os.system(plotcmd)

    # if only one directory, plot histogram without matches.
    if nDir == 1:
        nall = 0
        plotHistogram( directoryDate, calendar, nDir, rs, nDay, mjdRef, \
                       transitObject, raTransit, decTransit, \
                       eventDirs, nall, [ 0., 0.], [ 0, 0], [ 0., 0.], \
                       [0., 0.], [" "], 0)
        sys.exit()

    if doLog:
        logFile.close()
    print( "Count of Events and Event Groups: %d" % (nGroup))

    # now have a list of single events and multiple events, plot histogram
    
    plotHistogram( directoryDate, calendar, nDir, rs, nDay, mjdRef,
                   transitObject, raTransit, decTransit, \
                   eventDirs, nGroup, matchtimes, matchcounts, \
                   matchgallon, matchgallat, groupFlags, nGroup)

# now count all events happening within .1 degrees of other events.
    eventDaySummary( directoryDate, calendar, directoryDate + "/index.html")
    
if __name__ == "__main__":
    main()

