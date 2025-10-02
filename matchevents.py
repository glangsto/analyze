#Python find matchs in data directories
#HISTORY
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
import time
import datetime
import glob
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.coordinates import AltAz
import radioastronomy
import subprocess
import copy
from jdutil import date_to_mjd
from findTime import findTime
from findDate import findDirs, readEventsInDir
from findMatches import *
from groupMatches import *

nargs = len(sys.argv)
if nargs < 2:
    print("MATCH: Match events listed in several directories")
    print("Usage: MATCH [-OF seconds] [-D] [-C Date] [-E] [-N <n>] [-G] [dir1 dir2 dir3 dir4 ...]")
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
doQuiet = False
doHistogram = False
doLog = False
logFileName = ""
outdir = ""
flagGroups = False
doOffset = False
yoffset = 0.0
doGalactic = True
doGaussian = False
doTransit = False
raTransit = 0.0
decTransit = 0.0
nPrint = 4
minEl = -100.
# temporarily turn on/off debugging
verbose = True
verbose = False

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
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
    if str(anarg[0:3]) == "-W":
        outdir = sys.argv[iii+1]
        iii = iii + 1
        print("Writing plots to %s" % (outdir))
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
        raTransit = float(sys.argv[iii])
        iii = iii + 1
        decTransit = float(sys.argv[iii])
        print("Marking Transit of RA, Dec: %7.2f,%7.2f (deg)" % \
              (raTransit, decTransit))
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
        doOffset = True
        iii = iii + 1
        yoffset = float(sys.argv[iii])
        ifile = ifile + 2
    iii = iii + 1

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

# prepare to write out summary of events
if doLog:
    if os.path.exists(logFileName):
        # if an existing file
        if os.path.isfile(logFileName):
            logFile = open(logFileName, 'a')
        else:
            outName = logFileName + "/" + calendar + ".log"
            logFile = open(outName, 'w')
            print( "Writing to Log file: %s" % (outName))
    else:   # else neither path nor file exist, write a new file
        logFile = open(logFileName, 'w')
    print( "# Log of events on %s" % (calendar), file=logFile)
    print( "#  Tel   Directory             Event File            Ra     Dec     GLon    GLat N Flashes", file=logFile)
         
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
                

def plotHistogram( nDir, rs_in, nDay, mjdRef, doTransit, raTransit, decTransit, \
                   EventDirs, nall, matchtimes, matchcount, matchgallon, \
                   matchgallat, flags, nUnique):
    """
    plot several histograms of the event count versus time of day
    where:
    nDir = number of telescopes
    nUnique - number of distinque events
    """
    import subprocess
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    global calendar
    fig, ax = plt.subplots( figsize=(12,6))

    rs = copy.deepcopy(rs_in)

    # get local timezon and utc offset
    batcmd="/bin/date +%Z"
    timezone = subprocess.check_output(batcmd, shell=True)
    parts = timezone.split()
    # fix strange problem with zone sometimes coming back in quotes like b'EST'
    try:
        timezone = str(parts[0], 'UTF-8')
    except:
        timezone = str(parts[0])

    # prepart to compute local time
    utcOffset = datetime.datetime.utcnow() - datetime.datetime.now()
    utcOffsetSecs = utcOffset.total_seconds() + 1.
    utcOffsetSecs = int(utcOffsetSecs)
    utcOffsetHours = utcOffsetSecs/3600.
    utcparts = np.zeros(nDay+1)
    deltat = 24./float(nDay)
    for iDay in range(nDay):
        utcparts[iDay] = float(iDay*24./float(nDay))

    # to finish plot need to duplicate last point
    utcparts[nDay] = 24.
    if calendar == "":
        utcstr = str(rs.utc)
        strparts = utcstr.split()
        calendar = strparts[0]
        print( 'Found Calendar date: %s' % (calendar))
        
    countmax = np.zeros( nDir)
    yoffsets = np.zeros( nDir)
    yoffset = 0
    for ddd in range( nDir):
        counts = copy.deepcopy( EventDirs[ddd]['counts'])
        countmax[ddd] = np.max(counts)
        yoffset = yoffset + countmax[ddd]

    ytop = yoffset
    # now back down from the top for each plot
    for ddd in range(nDir):
        yoffset = yoffset - countmax[ddd]
        yoffsets[ddd] = yoffset

    yoffsets[nDir-1] = 0

#    barcolors = ["lightcyan", "wheat", "greenyellow", "mistyrose", \
#                 "wheat", "lightgrey", "mintcream"]
    barcolors = ["gold", "wheat", "greenyellow", "mistyrose", \
                 "wheat", "lightgrey", "mintcream"]
    # now for all directorys with events
    for iDir in range( nDir):
        adir = EventDirs[iDir]['dir']
        nEve = EventDirs[iDir]['n']
        counts = copy.deepcopy( EventDirs[iDir]['counts'])
#        if verbose:
#            for iDay in range(nDay):
#                if counts[iDay] > 0:
#                    print("%5d %8.3f %5d" % ( iDay, utcparts[iDay], counts[iDay]))

        alabel = adir[0:5]
        labelparts = alabel.split("-")
        alabel = "%4s:%4d" % (labelparts[0], nEve)
        # duplicate the last value to complete the plot
        counts = np.append( counts, counts[nDay-1])
        ybottom = yoffsets[iDir]
# time epsilon to unhide multiple events
        plt.step( utcparts, counts + yoffsets[iDir], where="post", label=alabel)
        plt.bar( utcparts, counts, width=deltat, bottom=ybottom, align="edge", \
                  label=None, color=barcolors[iDir])
        #                  where='post', label=alabel)
        maxcounts = np.max(counts)

# Redefine yoffset as step between printing lines
    yoffset = ytop / 30
    y0 = yoffset

# now annoated Az, El, Ra Dec
    alabel = " AZ,El:  %5.1f,%5.1f" % (rs.telaz, rs.telel)
    xa = utcOffsetHours
    ya = ytop
    ax.annotate(alabel, xy=( xa, ya), xytext=(xa , ya) )
    verticalcolors = ['pink','g','lightblue','c','m', 'k']
    # get Local Ra, Dec at UTC midnight
    utcstr = str(rs.utc)
    parts = utcstr.split()
    # get utc of begining of utc day
    dateHour = "%s %5d" % (parts[0], 0)
    utc0 = datetime.datetime.strptime( dateHour, "%Y-%m-%d %H")
    # get utc of local midnight
    dateHour = "%s %5d" % (parts[0],utcOffsetHours)
#    print("dateHour: %s" % (dateHour))
    utcmidnight = datetime.datetime.strptime( dateHour, "%Y-%m-%d %H")
    print("Utc: %s, Local Midnight Utc: %s" % (utcstr, utcmidnight))
    # prepart to compute local time
#    utcOffsetHours = time.timezone/3600.
    print( "Time zone %s is offset %5.1f hours from UTC" % (timezone, utcOffsetHours))

    # annoate plot for local midnight
    rs.utc = utcmidnight
    rs.azel2radec()
    xa = utcOffsetHours % 24.
    ramidnight = rs.ra   # keep ra (degrees) of midnight
    xmidnight = xa
    ya = ytop - yoffset
    alabel = " RA,Dec: %5.1f,%5.1f" % (rs.ra, rs.dec)
    ax.annotate(alabel, xy=( xa, ya), xytext=( xa, ya) )
    epsilon = 3./1440
    # now draw vertical lines for midnight local time
    plt.axvline( xa, color='b', linestyle='dotted')
    alabel = " Local Midnight"
    xa = utcOffsetHours % 24.
#    ybottom = -yoffset*.75
#    ax.annotate(alabel, xy=( xa, ybottom), xytext =( xa, ybottom))
    lst = datetime.timedelta(seconds=(86400.*rs.lst/360.))
    lststr = str(lst)
    lstparts = lststr.split(" ")
    nparts = len(lstparts)
    if nparts > 1:
        lststr = lstparts[nparts - 1]
    lstparts = lststr.split(".")
    lststr = lstparts[0]
# had trouble with LST offset across daylight savings time changes
#    alabel = " LST: %s" % (lststr)
#    ax.annotate(alabel, xy=( xa, yoffset*.25), xytext =( xa, yoffset*.25))

    # utc time of noon
    dateHour = "%s %5d" % (parts[0],utcOffsetHours+12.)
#    print("dateHour: %s" % (dateHour))
    utcnoon = datetime.datetime.strptime( dateHour, "%Y-%m-%d %H")
    rs.utc = utcnoon
    print( "          Utc Noon : %s" % (str(rs.utc)))
    # save telescope pointing direction
    telel = rs.telel
    telaz = rs.telaz
    # compute az, alt of sun at noon
    if telaz < 90. or telaz > 270.:
        rs.telaz = 180.
    # compute location of sun at noon
        rs.azel2radec()
    # then use el of sun at noon to compute ra,dec of sun
    rs.telaz = rs.az_sun
    rs.telel = rs.altsun
    rs.azel2radec()
    print( "Noon Sun Az, El    : %7.2f, %7.2f" % (rs.az_sun, rs.altsun))
    # numerically find sun transit time via newton's method
    for kkk in range(10):
        daz = (180. - rs.az_sun)/2.
#        print( "Daz %8.2f (rs.lst = %8.2f, rs.az_sun = %8.2f" % \
#           (daz, rs.lst, rs.az_sun))
        xtransit = 12. + (24*daz/360.) + utcOffsetHours
        rs.utc = rs.utc + datetime.timedelta(seconds=(86400.*daz/360.))
        rs.azel2radec()
    # Compute offset between noon sun Az and horn az (degrees)
    print( "Transit Sun        : %s" % (str(rs.utc)))
    print( "Transit Sun Az, El : %7.2f, %7.2f" % (rs.az_sun, rs.altsun))
    rs.telaz = rs.az_sun
    rs.telel = rs.altsun
    rs.azel2radec()
    transit = rs.utc - utc0
#    print("dt: %s" % transit)
    dt = transit.total_seconds()/86400.
    xtransit = 24.*dt
#    print("dt: %s dt=%7.3f  %7.3f" % (transit, dt, xtransit))
    alabel = " RA,Dec: %5.1f, %5.1f" % (rs.ra, rs.dec)
    ya = ytop
    xa = xtransit
    ax.annotate(alabel, xy=( xa, ya), xytext =( xa, ya))
    # set elevation back to horn elevaiton
    alabel = " AZ,El: %5.1f, %5.1f" % (rs.az_sun, rs.altsun)
    ya = ytop - yoffset
    ax.annotate(alabel, xy=( xtransit, ya), xytext =( xtransit, ya))
    # restore telescope az and el
    rs.telaz = telaz
    rs.telel = telel
    alabel = "Sun:"
    ya = ytop - (0.5*yoffset)
    xa = xtransit - 1.
    ax.annotate(alabel, xy=( xa, ya), xytext =( xa, ya))

    # now draw vertical lines for Sun cross horn
    plt.axvline( xtransit, color='b', linestyle='dotted')


    # next draw galactic az - 10 and galactic z + 10 vertial lines
    lat = rs.tellat
    lon = rs.tellon
    try:
        height = rs.telelev
    except:
        height = 1000.
    # first need site location to get coordinates of galaxy
    asite = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    tutc = utc0   # start at beginning of day to find transit utc
    target = SkyCoord( ra=raTransit*u.deg, dec=decTransit*u.deg,
                       unit='deg', \
                       frame='icrs', location=asite, obstime=tutc)
    # must convert ras from degrees to hours
    xtransit = ((raTransit - ramidnight)/15.) + xmidnight
    if doTransit:
        tutc = utc0   # start at beginning of day to find transit utc
        dutc = 0
        if xtransit > 24.:
            xtransit = xtransit - 24.
        if xtransit > 24.:
            xtransit = xtransit - 24.
        if xtransit < 0.:
            xtransit = xtransit + 24.
        if xtransit < 0.:
            xtransit = xtransit + 24.
        # now draw vertical lines for Sun cross horn
        print( "X transit = %8.2f" % (xtransit))
        plt.axvline( xtransit, color='r', linestyle='dashed')
        alabel = " RA,Dec: %5.1f, %5.1f" % (raTransit, decTransit)
        ya = ytop - yoffset - yoffset
        ax.annotate(" Transit ", xy=( xtransit, ya), xytext =( xtransit, ya))
        alabel = " RA,Dec: %5.1f, %5.1f" % (raTransit, decTransit)
        ya = ya - yoffset
        ax.annotate(alabel, xy=( xtransit, ya), xytext =( xtransit, ya))

    # prepare to draw +- 10 degrees of galactic plane
    nGalactic = 90
    nG4 = int(nGalactic/4)
    if doGalactic:
        utc = utc0
        dt = datetime.timedelta(seconds=(86400./nGalactic))
        minlat = 10.
        for mmm in range(nGalactic):   # every so often during of the day
            azel = SkyCoord(az = float(rs.telaz)*u.deg, \
                            alt = float(rs.telel)*u.deg, unit='deg', \
                            frame='altaz', location=asite, obstime=utc)
            if mmm % nG4 == 0:
                print( "MMM: %4d  %7.2f, %7.2f; %s" % \
                           (mmm, azel.galactic.l.degree, azel.galactic.b.degree, utc))
            bbb = azel.galactic.b.degree
            lll = azel.galactic.l.degree
            utc = utc + dt
            ygal = yoffset/2.    # set gold galactic line placement at top
            if bbb < minlat and bbb > -minlat:
                xgal = 24.*mmm/nGalactic
                ax.annotate("*", xy=( xgal, ygal), color='orange')
                ax.annotate("*", xy=( xgal, yoffset), color='orange')
                ax.annotate("*", xy=( xgal, ya-(1.5*yoffset)), color='orange')
                ax.annotate("*", xy=( xgal, ya-(2.0*yoffset)), color='orange')

    # now draw vertical lines for events
    iplot = 0
    # for all matched events
    
    for iii in range(nall):
        if flagGroups:    # If not flagging groups, skip this processing
            continue
        x4 = matchtimes[iii] - mjdRef
        # if date is beyond 0-24 hour range, move in range
        x4 = (x4*24.) % 24.
        # now compute x position for this MJD
        # MJD midnight = UTC midnight.
        # For EST, the time is actually 5 hours before midnight
        x4 = x4 + (epsilon*((iii%7)-3))
        icolor = iplot % 6
        plt.axvline( x4, color=verticalcolors[icolor], linestyle='dashed')
        # if a single match
        if flags[iii] == " ":
            alabel = "o %.1f,%.1f" % ( matchgallon[iii], matchgallat[iii])
        elif  flags[iii] == "X" or flags[iii] == "C" or flags[iii] == "M":
            alabel = "%s %.1f,%.1f" % (flags[iii], matchgallon[iii], matchgallat[iii])
        else:
            alabel = "%sx %.1f,%.1f" % (flags[iii], matchgallon[iii], matchgallat[iii])

        if flags[iii] == "M":
            ax.annotate(alabel, xy=( x4, y0), xytext=(x4 , y0), fontsize=12, color='red')
        elif flags[iii] == "C":
            ax.annotate(alabel, xy=( x4, y0), xytext=(x4 , y0), fontsize=11, color='green')
        elif flags[iii] == "X":
            ax.annotate(alabel, xy=( x4, y0), xytext=(x4 , y0), fontsize=10, color='blue')
        else:
            ax.annotate(alabel, xy=( x4, y0), xytext=(x4 , y0), fontsize=9, color='black')
        y0 = y0 + yoffset
        iplot = iplot + 1
        # keep from going off the plot
        if y0 > ytop - (2 * yoffset):
            y0 = 2 * yoffset
    # now finish up labels
    plt.legend()
    if nDir == 1:
        plt.title("%s   Event Rate" % (calendar),fontsize=18)
        plt.xlabel("Time (UTC hours)", fontsize=18)
    else:
        plt.title("%s   Event Rate for %d Observations" % (calendar, nDir),
                  fontsize=18)
#        plt.xlabel("Time (UTC hours, %s offset: %2.0f hours) (Distinct Event Count: %d)" % (timezone, utcOffsetHours, nUnique), fontsize=18)
        plt.xlabel("Time (UTC hours) (Distinct Event Count: %d)" % (nUnique), fontsize=18)

    if nDay == 24:
        plt.ylabel("Count of Events/Hour", fontsize=18)
    elif nDay == 1440:
        plt.ylabel("Count of Events/Minute", fontsize=18)
    else:
        nminute = int(1440.001/float(nDay))
        plt.ylabel("Count of Events/%d Minutes" % (nminute), fontsize=18)

    plt.xlim(0.,24.05) # 24 hours/per day
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    xticks = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
#    plt.xticks( float( np.array(xticks)), str(np.array(xticks)))
    nticks = len(xticks)
    strticks = xticks
    for i in range(nticks):
        strticks[i] = str( xticks[i])
        xticks[i] = float(xticks[i])
    ax.set_xticks( xticks)
   
    #    ax.set_xticklabels( strticks)
    ax.tick_params(which='major', width=3)
    ax.tick_params(which='minor', width=1)
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=3, color='b')
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.25))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    if outdir != "":
        plotfile = '%s/match-%s.pdf' % (outdir, calendar)
    else:
        plotfile = 'match-%s.pdf' % (calendar)
    fig = plt.gcf()
    fig.savefig(plotfile, bbox_inches='tight')
    print("Saved match summary plot: %s" % (plotfile))
    if outdir != "":
        plotfile = '%s/match-%s.svg' % (outdir, calendar)
    else:
        plotfile = 'match-%s.svg' % (calendar)
    fig = plt.gcf()
    fig.savefig(plotfile, bbox_inches='tight')
    if not doQuiet:
        plt.show()
    return
# end of plotHistogram

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
        if mjdRef == 0:
            mjdRef = int(mjds[0])
        if int(mjds[0]) < mjdRef:
            mjdRef = int(mjds[iDir])

    if nDir > 0:
        print( "MJD reference: %12.6f" % (mjdRef))

    rs = radioastronomy.Spectrum()
    nTotal = 0   # count total number of events in all directories
    # next, for each directory, count eventNames in each hour/part of a day
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
        

    matchCounts = np.zeros( nDir)

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
    nGroup, nTen, nHundred, nThousand, groupEvents = groupMatches( mjdRef, \
                                    nRemain, nDir, eventTrim, nDay, verbose = True)
    print("Grouped %d events into %d Groups" % \
          (nRemain, nGroup))
    print("In Groups, %3d with > 10 members, %d with > 100 and %d > 1000 members" % \
          (nTen, nHundred, nThousand))

    rs = radioastronomy.Spectrum()

    matchtimes = np.zeros(nGroup)
    matchcounts = np.zeros(nGroup)
    matchgallon = np.zeros(nGroup)
    matchgallat = np.zeros(nGroup)

#    verbose = False

    verbose = True
    
    # initialize the file list
    files = []
    groupFlags = [ ]
    # now fill arrays with coordinates and info
    for iGroup in range(nGroup): 
        aGroup = groupEvents[ iGroup] 
        groupFlags.append( aGroup['flag'])
        matches = aGroup['list']
        # count becomes number of events in a group
        nEvents = int(aGroup['count'])
        matchtimes[iGroup] = aGroup['mjd']
        matchcounts[iGroup] = nEvents
        files = []
        nTel = 0
        aveGalLon = 0.
        aveGalLat = 0.
        for iDir in range(nDir):
            iMatch = matches[iDir]
            if iMatch != NOMATCH:
                # find file name in sets of names
                fileName = eventDirs[iDir]['events'][iMatch]
                if doPlot:
                    files.append( fileName)
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
            
        if verbose:
            aveMjd = aGroup['mjd']
            showMatch( iGroup, aveMjd, nDir, matches, eventDirs)
                             
        if doPlot:
            if not doOffset:
                yoffset = 0.
            
            plotcmd = "~/Research/analyze/E -Y %.2f %s " % (yoffset, files)
            os.system(plotcmd)
            
               
    # if only one directory, plot histogram without matches.
    if nDir == 1:
        plotHistogram( nDir, rs, nDay, mjdRef, doTransit, \
                       raTransit, decTransit, \
                       eventDirs, 0, [ 0., 0.], [ 0, 0], [ 0., 0.], [0., 0.], [" "], 0)
        exit()

    if doLog:
        logFile.close()
    print( "Count of Events and Event Groups: %d" % (nGroup))

# now have a list of single events and multiple events, plot histogram
    plotHistogram( nDir, rs, nDay, mjdRef, doTransit, \
                   raTransit, decTransit, \
                   eventDirs, nGroup, matchtimes, matchcounts, \
                   matchgallon, matchgallat, groupFlags, nGroup)

# now count all events happening within .1 degrees of other events.

if __name__ == "__main__":
    main()

