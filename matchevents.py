#Python find matchs in data directories
#HISTORY
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

nargs = len(sys.argv)
if nargs < 2:
    print("MATCH: Match events listed in several directories")
    print("Usage: MATCH [-OF seconds] [-D] [-C Date] [-E] [-N <n>] [-G] [dir1 dir2 dir3 dir4 ...]")
    print("Where: Optionally the user provides the maximum time offset (secs) to call a match")
    print("Where -C  Optionally provide a calendar date (ie 19Nov17) instead of directories")
    print("Where -D  Optionally print debugging info")
    print("Where -E Optionally only show Events with elevation above zero")
    print("Where -G Optionally plot +/- 10 degrees of galactic plane")
    print("Where -H Optionally plot histogram of events per parts of a day")
    print("Where -N <n> Optionally print matches when number is equal or greater to <n>")
    print("Where -P Plot matches")
    print("Where -Y <offset> Optionally adds an offset to the event plots (-P option)")
    print("")
    print("Glen Langston, 2021 September 1")
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
doOffset = False
yoffset = 0.0
doGalactic = False
nPrint = 4
minEl = -100.

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-OF":
        offset = float(sys.argv[iii+1])
        iii = iii + 1
        print("Maximum Time Offset: %8.6f s for Match: " % (offset))
        offset = offset/86400.   # convert to MJDs
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
        print("Only counting events when telescope above %5.1f (d) elevation" % (minEl))
        ifile = ifile + 1
    if str(anarg[0:3]) == "-F":
        flagGroups = True
        print("Flagging groups of events")
        ifile = ifile + 1
    if str(anarg[0:3]) == "-G":
        doGalactic = True
        print("Marking observations within 10 degrees of the galatic plane with *")
        ifile = ifile + 1
    if str(anarg[0:3]) == "-ND":
        nday = int(sys.argv[iii+1])
        iii = iii + 1
        print("Divide Day into N Parts:  %d" % (nday))
        aFix = True
        ifile = ifile + 2
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
    if anarg[0:2] == "-H":
        doHistogram = True
        print("Plotting Histogram of events per day")
        ifile = ifile + 1
    if anarg[0:2] == "-D":
        doDebug = True
        print("Debugging")
        ifile = ifile + 1
    if anarg[0:2] == "-Y":
        doOffset = True
        iii = iii + 1
        yoffset = float(sys.argv[iii])
        ifile = ifile + 2
    iii = iii + 1

nplot = 0
if nday < 1:
    nday = 24

#print("Counting Events in blocks of %5.2f hours" % (24./nday))
# 
nDir = nargs-ifile
# 
MAXDIR = 10
NOMATCH = -9999
MAXDT = 10000.
# define the maximum number of events to keep
MAXEVENTS = 10000


def finddirs( calendar, ifile, nDir):
    """
    finddirs() finds the names of directories matching the calendar flag
    """
    dirs = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""]
    nDir = 0
    nfiles = len(sys.argv[ifile:])
    # now find calendars
    if calendar == "":
        dirs[0] = sys.argv[ifile]
        if nfiles < 2:
            nDir = 1
            return dirs[0:1]
        dirs[1] = sys.argv[ifile+1]
        nDir = 2
        if n > 2:
            dir[nDir] = sys.argv[ifile+2]
            nDir = nDir+1
        if nDir > 3:
            dir[nDir] = sys.argv[ifile+3]
            nDir = nDir+1
    else:
        for idir in range(1, 40):
            adir = "pi%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs[nDir] = adir
                nDir = nDir + 1
                print("Found directory: %s" % adir)
            if nDir >= MAXDIR:
                break

        for idir in range(1, 40):
            adir = "odroid%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs[nDir] = adir
                nDir = nDir + 1
                print("Found directory: %s" % adir)
            if nDir >= MAXDIR:
                break

        for idir in range(1, 20):
            adir = "od%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs[nDir] = adir
                nDir = nDir + 1
                print("Found directory: %s" % adir)
            if nDir >= MAXDIR:
                break

    if nDir < 2:
        print("Can not match events in less than 2 directories")

    return dirs[0:nDir]

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
        if rs.telel < minEl:
            continue
        mjds[kkk] = rs.emjd
        peaks[kkk] = rs.epeak
        rmss[kkk] = rs.erms
        azs[kkk] = rs.telaz
        els[kkk] = rs.telel
        if (100 * int(kkk/100) == kkk) and doDebug:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % (kkk, rs.emjd, rs.epeak, rs.erms))
        kkk = kkk + 1
    nEve = kkk

    return nEve, events, mjds, peaks, rmss, azs, els
# end of readEventsInDir

# count directories with events
dirs = finddirs( calendar, ifile, nDir)

if nday != 24:
    print("Counting Events in blocks of %5.2f hours" % (24./nday))

if doDebug:
    nDir = len(dirs)
    for i in range(ndir):
        print("Dir %d: %s" % (i+1, dirs[i]))

def findpairs( EventDir1, EventDir2):
    """
    findpairs() - match events in pairs of event structures
    Inputs are Python data structures describing the events in the directories
    Outputs are indicies of matches to 2nd directory and time offsets between
    """
    # extract number of events in each directory
    n1 = EventDir1['n']
    n2 = EventDir2['n']
    # make sure there are plenty of indicies, set to no-match
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

def plotHistogram( nDir, rs_in, nday, mjdRef, EventDirs, nall, match4times, match4count, match4gallon, match4gallat, nUnique):
    """
    plot several histograms of the event count versus time of day
    """
    import subprocess
    import matplotlib.pyplot as plt
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
    utcparts = np.zeros(nday+1)
    for iDay in range(nday):
        utcparts[iDay] = float(iDay*24./float(nday))
    # to finish plot need to duplicate last point
    utcparts[nday] = 24.
    if calendar == "":
        utcstr = str(rs.utc)
        strparts = utcstr.split()
        calendar = strparts[0]

    yoffset = 0
    barcolors = ["lightcyan", "wheat", "greenyellow", "mistyrose", \
                 "wheat", "lightgrey", "mintcream"]
    for ddd in range( nDir):
        adir = EventDirs[ddd]['dir']
        nEve = EventDirs[ddd]['n']
        counts = copy.deepcopy( EventDirs[ddd]['counts'])
        alabel = adir[0:5]
        labelparts = alabel.split("-")
        alabel = labelparts[0]
        # duplicate the last value to complete the plot
        counts = np.append( counts, counts[nday-1])
        ybottom = yoffset + (0.* counts)
# time epsilon to unhide multiple events
        plt.step( utcparts, counts + yoffset, where="post", label=alabel)
#        plt.bar( utcparts, counts, bottom=ybottom, align="center", \
#        plt.bar( utcparts, counts, bottom=ybottom, align="edge", \
        plt.bar( utcparts, counts, bottom=ybottom, align="edge", \
                  label=None, color=barcolors[ddd])
        #                  where='post', label=alabel)
        maxcounts = np.max(counts)
        yoffset = yoffset + maxcounts
# keep the maximum of all points for scaling plot labels
    max0 = yoffset
# put 20 annotation lines vertically
    yoffset = yoffset / 20
    y0 = 2*yoffset

# now annoated Az, El, Ra Dec
    alabel = " AZ,El:  %5.1f,%5.1f" % (rs.telaz, rs.telel)
    xa = utcOffsetHours
    ya = max0 - yoffset
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
    alabel = " RA,Dec: %5.1f,%5.1f" % (rs.ra, rs.dec)
    xa = utcOffsetHours % 24.
    xmidnight = xa
    ytop = max0
    ya = ytop
    ax.annotate(alabel, xy=( xa, ya), xytext=( xa, ya) )
    epsilon = 3./1440
    # now draw vertical lines for midnight local time
    plt.axvline( xa, color='b', linestyle='dotted')
    alabel = " Local Midnight"
    xa = utcOffsetHours % 24.
    ybottom = -yoffset*.75
    ax.annotate(alabel, xy=( xa, ybottom), xytext =( xa, ybottom))
    lst = datetime.timedelta(seconds=(86400.*rs.lst/360.))
    lststr = str(lst)
    lstparts = lststr.split(" ")
    nparts = len(lstparts)
    if nparts > 1:
        lststr = lstparts[nparts - 1]
    lstparts = lststr.split(".")
    lststr = lstparts[0]
    alabel = " LST: %s" % (lststr)
    ax.annotate(alabel, xy=( xa, yoffset*.25), xytext =( xa, yoffset*.25))

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
    ya = ytop-yoffset
    ax.annotate(alabel, xy=( xtransit, ya), xytext =( xtransit, ya))
    # restore telescope az and el
    rs.telaz = telaz
    rs.telel = telel
    ya = ytop
    xa = xtransit
    alabel = "Sun:"
    ya = ytop - (yoffset/2.)
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
    # first need site clcation to get coordinates of galaxy
    asite = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)

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
            ygal = ytop * 1.02    # set gold galactic line placement at top
            if bbb < minlat and bbb > -minlat:
                xgal = 24.*mmm/nGalactic
                ax.annotate("*", xy=( xgal, ygal), color='gold')

    # now draw vertical lines for events
    iplot = 0
    for iii in range(nall):
        x4 = match4times[iii] - mjdRef + (0.5/float(nday))
#        x4 = match4times[iii] - mjdRef
        x4 = (x4*24.) % 24.
        # now compute x position for this MJD 
        # MJD midnight = UTC midnight.
        # For EST, the time is actually 5 hours before midnight
        x4 = x4 + (epsilon*((iii%7)-3))
        icolor = iplot % 6
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
        if match4count[iii] >= 1:
            y0 = y0 + yoffset
            iplot = iplot + 1
            # keep from going off the plot
        if y0 > max0 - (2 * yoffset):
            y0 = 2 * yoffset
    # now finish up labels
#    plt.legend(title="Set of Obs.")
    plt.legend()
    if nDir == 1:
        plt.title("%s   Events per Hour" % (calendar))
        plt.xlabel(("Time (UTC hours, %s offset: %2.0f hours)" % (timezone, utcOffsetHours)), fontsize=16)
    else:
        plt.title("%s   Events per Hour for %d Observations" % (calendar, nDir))
        plt.xlabel("Time (UTC hours, %s offset: %2.0f hours) (Distinct Event Count: %d)" % (timezone, utcOffsetHours, nUnique), fontsize=16)
    plt.ylabel("Count of Events/Hour", fontsize=16)
    plt.xlim(0.,24.) # 24 hours/per day
    plotfile = 'match-%s.png' % (calendar)
    fig = plt.gcf()
    fig.savefig(plotfile, bbox_inches='tight')
    
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

    nDir = len(dirs)
    nEvents = np.zeros(nDir)
          
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

    print( "MJD reference: %12.6f" % (mjdRef))

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
            dMjd = mjds[iEve] - mjdRef + (0.5/float(nday))
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
    if nDir > 5:
        session5 = EventDirs[ 5 ]
        nEve5 = EventDirs[ 5 ]['n']
        counts5 = copy.deepcopy( EventDirs[5]['counts'])
    else:
        nEve5 = 0
        counts5 = np.zeros(nday)

    print("Hour        Telescope/Day  ")
    if nDir > 4:
        print("Hour      1     2     3     4     5")
    elif nDir > 3:
        print("Hour      1     2     3     4")
    elif nDir > 2:
        print("Hour      1     2     3")
    else:
        print("Hour      1     2")

    utcparts = np.zeros(nday)
    # now print matches as a function of time 
    for iDay in range(nday):
        utcparts[iDay] = float(iDay*24./float(nday))
        if counts0[iDay] != 0 or counts1[iDay] != 0 or counts2[iDay] != 0 \
                or counts3[iDay] != 0 or counts4[iDay] != 0:
            if nDir == 5:
                print("%5.1f  %5d %5d %5d %5d %5d" % (utcparts[iDay], \
                                                  counts0[iDay],counts1[iDay], counts2[iDay], counts3[iDay], counts4[iDay]))
            elif nDir == 4:
                print("%5.1f  %5d %5d %5d %5d" % (utcparts[iDay], \
                                               counts0[iDay],counts1[iDay], counts2[iDay], counts3[iDay]))
            elif nDir == 3:
                print("%5.1f %5d %5d %5d " % (utcparts[iDay], \
                                               counts0[iDay],counts1[iDay], counts2[iDay]))
            elif nDir == 2:
                print("%5.1f %5d %5d " % (utcparts[iDay], \
                                               counts0[iDay],counts1[iDay]))
            elif nDir == 1:
                print("%5.1f %5d " % (utcparts[iDay], \
                                               counts0[iDay]))
            else: # else 5 or more, just print 5
                print("%5.1f %5d %5d %5d %5d %5d" % (utcparts[iDay], counts0[iDay], \
                                                       counts1[iDay], counts2[iDay], counts3[iDay], counts4[iDay]))

# Given N telescopes, there are N!/2 pairs of observations 
# ie 2 Telescopes, 1 Pair
#    3 Telescopes, 3 Pairs
#    4 Telescopes  6 Pairs
#    5 Telescopes 30 Pairs

# determine the maximum number of matches
    maxMatch = nTotal 

    mjd0s = EventDirs[0]['mjds']
    # for each of the pairs of directories, find closest matches
    if nDir > 1:
        mjd1s = EventDirs[1]['mjds']
        ii01s, dt01s = findpairs( EventDirs[0], EventDirs[1])
    if nDir > 2:
        mjd2s = EventDirs[2]['mjds']
        ii02s, dt02s = findpairs( EventDirs[0], EventDirs[2])
        ii12s, dt12s = findpairs( EventDirs[1], EventDirs[2])
    if nDir > 3:
        mjd3s = EventDirs[3]['mjds']
        ii03s, dt03s = findpairs( EventDirs[0], EventDirs[3])
        ii13s, dt13s = findpairs( EventDirs[1], EventDirs[3])
        ii23s, dt23s = findpairs( EventDirs[2], EventDirs[3])
    if nDir > 4:
        mjd4s = EventDirs[4]['mjds']
        ii04s, dt04s = findpairs( EventDirs[0], EventDirs[4])
        ii14s, dt14s = findpairs( EventDirs[1], EventDirs[4])
        ii24s, dt24s = findpairs( EventDirs[2], EventDirs[4])
        ii34s, dt34s = findpairs( EventDirs[3], EventDirs[4])

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
        i4 = NOMATCH
        if i0 == 0:
            match = { 'nmatch': nMatch, 'count': matchcount, 'mjd': mjdRef, 'list': [ i0, i1, i2, i3] }
            matchs = { nMatch: match }
            
            
        if nDir < 2:
            continue
        if abs(dt01s[i0]) < offset:
            i1 = int(ii01s[i0])
            mjdave = mjdave + mjd1s[i1]
            matchcount = matchcount + 1
        if nDir > 2:
            if abs(dt02s[i0]) < offset:
                i2 = int(ii02s[i0])
            #            print("0: %5d, 1: %5d 2: %5d" % (i0, i1, i2))
                mjdave = mjdave + mjd2s[i2]
                matchcount = matchcount + 1
        if nDir > 3:
            if abs(dt03s[i0]) < offset:
                i3 = int(ii03s[i0])
                mjdave = mjdave + mjd3s[i3]
                matchcount = matchcount + 1
        if nDir > 4:
            if abs(dt04s[i0]) < offset:
                i4 = int(ii04s[i0])
                mjdave = mjdave + mjd4s[i4]
                matchcount = matchcount + 1
        # if any matches to this event
        if matchcount > 1:
            # a match has at least two partners
            mjdave = mjdave / float(matchcount)
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
        if nDir < 2:
            continue
        if abs(dt12s[i1]) < offset:
            i2 = int(ii12s[i1])
            mjdave = mjdave + mjd2s[i2]
            matchcount = matchcount + 1
        if nDir > 3:
            if abs(dt13s[i1]) < offset:
                i3 = int(ii13s[i1])
                mjdave = mjdave + mjd3s[i3]
                matchcount = matchcount + 1
        if matchcount > 1:
            match['nmatch'] = nMatch
            match['count'] = matchcount
            match['mjd'] = mjdave / float(matchcount)
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
        if nDir > 3:
            if abs(dt23s[i2]) < offset:
                i3 = int(ii23s[i2])
                mjdave = mjdave + mjd3s[i3]
                matchcount = matchcount + 1
        if matchcount > 1:
            match['nmatch'] = nMatch
            match['count'] = matchcount
            match['mjd'] = mjdave / float(matchcount)
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
        if nDir < 2:
            continue
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
        if nDir < 2:
            continue
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
    print( " Tel:     2     3     4     5")
    nCheck = min( nDir, 4)
    for iii in range(nCheck): 
        counttypes = np.zeros(max(nDir,4))
        for lll in range(nMatch):
            matcha = matchs[ lll]
            lista = matcha['list']
            counta = matcha['count'] - 1
            if counta < 0:
                continue
            if lista[iii] != NOMATCH:
                counttypes[counta] = counttypes[counta] + 1
        print( "%3d:  %5d %5d %5d" % (iii, counttypes[1], counttypes[2], counttypes[3]))

    # prepare a full spectrum for time reference
    rs = radioastronomy.Spectrum()
    # read in the first event
    fullname = EventDirs[ 0]['events'][0]
    rs.read_spec_ast( fullname)
    rs.azel2radec()

    if nDir == 1:
        plotHistogram( nDir, rs, nday, mjdRef, EventDirs, 0, [ 0., 0.], [ 0, 0], [ 0., 0.], [0., 0.], 0)
        exit() 

#    for iii in range(nDir): 
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
#        print( "Directory: %2d - %5d"  % (iii, counts[iii]))

    nall=0
    # make big array to keep track of multiple matches
    files0 = EventDirs[0]['events']
    match4times = np.zeros(nMatch)
    match4gallon = np.zeros(nMatch)
    match4gallat = np.zeros(nMatch)
    match4index = np.zeros(nMatch)
    # now make a list of events that only have enough matches
    for lll in range(nMatch):
        matcha = matchs[ lll ]
        lista = matcha['list']
        counta = matcha['count']
        if counta < nPrint:
            continue
        # must have one directory, but others might not be present
        i0 = lista[0]
        if i0 < 0:
            continue
        if i0 >= nEve0:
            print("Error reading telescope 0 index %d not in range 0 to %d" % \
                  ( i0, len(files0)))
            print(lista)
            continue
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
    nUnique = 0
    # next print complete matches within the same (small) region
    for lll in range(nMatch):
        matcha = matchs[ lll ]
        lista = matcha['list']
        counta = matcha['count']
        if counta < nPrint:
            continue
        # now only print unique or 1st of many matches
        isUnique = False
        for iii in range(nall):
            if lll == match4index[iii]:
                isUnique = True
                multiples = match4count[iii]
        if not isUnique:
            continue
        nUnique = nUnique + 1
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
            rs.azel2radec()
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
            file4 = files4[i4]
            print("%3d 4 %s" %  (lll, file4))
        if doPlot:
            if doOffset:
                plotcmd = "~/Research/analyze/E -Y %.2f %s %s %s %s %s" % (yoffset, file0, file1, file2, file3, file4)
            else:
                plotcmd = "~/Research/analyze/E %s %s %s %s %s" % (file0, file1, file2, file3, file4)
            os.system(plotcmd)

    print( "Count of Events and Event Groups: %d" % (nUnique))

# now have a list of single events and multiple events, plot histogram
    plotHistogram( nDir, rs, nday, mjdRef, EventDirs, nall, match4times, match4count, match4gallon, match4gallat, nUnique)
    
# now count all events happening within .1 degrees of other events.
        
if __name__ == "__main__":
    main()

