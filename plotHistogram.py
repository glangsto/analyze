#Python find matchs in data directories and write html
#HISTORY
#25Dec15 GIL restore transit
#25Dec13 GIL separate out plotHistogram to s separate file
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
import utcOffset


# Only read one event file, to get location and az, el
readAzEl = False

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def plotHistogram( directoryDate, calendar, nDir, rs_in, nDay, mjdRef, \
                   transitObject, raTransit, decTransit, \
                   EventDirs, nall, matchtimes, matchcount, \
                   matchgallon, matchgallat, flags, nUnique):
    """
    plot several histograms of the event count versus time of day
    where:
    directoryDate = output directory for histograms
    calendar  = calendar date of events (ie 25Dec15)
    nDir      = number of telescopes
    rs_in     = radio astronomy structure containing one event
    nDay      = number of divisions of a day for histogram
    mjdRef    = reference mjd, often the first mjd of the events
    transitObject = name of object to mark RA, Dec of Transkit
    raTransit, decTransit = coordinate (degrees) of object, used to calc transit
    EventDirs = names of event directories (ie pi1-events-25Dec15)
    nall      = number of events to be marked
    matchtimes= times of events to be marked
    matchcount= Number of flashes seen for this event
    matchgallon, gallat = galactic coordinates of event
    flags     = single characters marking the number of flashes for this event
    nUnique   = number of distinct events
    
    """
    fig, ax = plt.subplots( figsize=(12,6))

    rs = copy.deepcopy(rs_in)
    utc = rs.utc
    utcStr = utc.strftime("%Y-%m-%d")
    # prepart to compute local time
    utcOffsetHours = - utcOffset.getUtcOffset( utcStr)

    utcparts = np.zeros(nDay+1)
    # now prepare to plot
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

    barcolors = ["gold", "wheat", "greenyellow", "mistyrose", \
                 "wheat", "lightgrey", "mintcream"]
    # now for all directorys with events
    for iDir in range( nDir):
        # adir - last part of file name before event. ie pi1-events-25Dec15
        adir = EventDirs[iDir]['dir']
        nEve = EventDirs[iDir]['n']
        counts = copy.deepcopy( EventDirs[iDir]['counts'])
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
#    raTransit = rs.ra
#    decTransit = rs.dec
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
    ya = ytop - yoffset
    xa = xtransit
    ax.annotate(alabel, xy=( xa, ya), xytext =( xa, ya))
    # set elevation back to horn elevaiton
    alabel = " AZ,El: %5.1f, %5.1f" % (rs.az_sun, rs.altsun)
    ya = ytop
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
    # must convert ras from degrees to hours
    xtransit = ((raTransit - ramidnight)/15.) + xmidnight
    if transitObject != "":
        tutc = utc0   # start at beginning of day to find transit utc
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
        ax.annotate(" Transit %s" % (transitObject), xy=( xtransit, ya), xytext =( xtransit, ya))
        alabel = " RA,Dec: %5.1f, %5.1f" % (raTransit, decTransit)
        ya = ya - yoffset 
        ax.annotate(alabel, xy=( xtransit, ya), xytext =( xtransit, ya))

    # prepare to draw +- 10 degrees of galactic plane
    nGalactic = 90
    nG4 = int(nGalactic/4)
    doGalactic = True
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
        plt.title("%s   Event Rate" % (calendar), fontsize=18)
        plt.xlabel("Time (UTC hours)", fontsize=18)
    else:
        plt.title("%s   Event Rate for %d Observations" % (calendar, nDir),
                  fontsize=18)
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
    plotfile = '%s/match-%s.pdf' % (directoryDate, calendar)
    fig = plt.gcf()
    fig.savefig(plotfile, bbox_inches='tight')
    print("Saved match summary plot: %s" % (plotfile))
    if directoryDate != "":
        plotfile = '%s/match-%s.svg' % (directoryDate, calendar)
    else:
        plotfile = 'match-%s.svg' % (calendar)
    fig = plt.gcf()
    fig.savefig(plotfile, bbox_inches='tight')
#    if not doQuiet:
#        plt.show()
    # return name of summary directory
    return 
# end of plotHistogram

# prepare to write out summa5ry of events
def logDirectory( calendar, outDirName, eventLogName):
    """
    logDirectory() opens files and creates paths so summarize events for a date
    where
    calendar   - calendar name of the events to summarize (25Dec12)
    outDirName - path to summary directory.  New subdirectories added here
    eventLogName - full path to single file loging all events file
    be created for month of year and day of year
    Returns pointer to open log file and and name of directory for
    writing one day of events.
    """
    from pathlib import Path

    expandedPath = Path(eventLogName).expanduser()
    pathDirectory = expandedPath.parent
    pathDirectory.mkdir(parents=True, exist_ok=True)
    logFileName = expandedPath.resolve()
    logFileName = str( logFileName)
    
    # if an existing file
    if os.path.isfile(logFileName):
        logFile = open(logFileName, 'a')  # append
    else:
        logFile = open(logFileName, 'w')  # write to a new file
    print( "Writing to Log file: %s" % (logFileName))

    expandedPath = Path(outDirName).expanduser()
    expandedPath.mkdir(parents=True, exist_ok=True)
    expandedOut = expandedPath.resolve()
    expandedOut = str( expandedOut)
    # now create month and day directories
    yearDir = "20" + calendar[0:2]
    monthDir = calendar[0:5]
    dateDir = calendar
    # create the full path
    directoryDate = expandedOut + "/" + yearDir + "/" + monthDir + "/" + dateDir
    os.makedirs(directoryDate, exist_ok=True)

    print("DirecotoryDate: %s" % (directoryDate))
    # return logFile pointer and new directory path
    return logFile, directoryDate

def main():
    """
    Test program for plotHistogram.
    The test requires an edit for you system, to define on event file
    location for reading site coordinates.
    """
    import subprocess

    # these parameters must be updated for the test to work on your system
    calendar = "25Dec13"
    HOME = "~/"
    dataDir = "/media/karl/"
    aDir = "pi1-events-25Dec13"
    anEvent = dataDir + aDir + '/25-12-13T153613_724.eve'
    eventDirs = ["pi1-events-25Dec13", \
                 "pi2-events-25Dec13", \
                 "pi3-events-25Dec13"]

    # end of your updates to test plotHistogram.

    # these directories and files will be created on your system
    outDirName = HOME + "testmatch"
    eventLogName = HOME + "testmatch/event.log"

    nDir = 3
    nGroup = 2
    matchtimes = np.zeros(nGroup)
    matchcounts = np.zeros(nGroup)
    matchgallon = np.zeros(nGroup)
    matchgallat = np.zeros(nGroup)

    
    # read in your event file here
    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(anEvent)
    rs.azel2radec()
    
    nDay = 144
    mjdRef = int(rs.emjd)
    mjds = [ rs.emjd]
    countsInDt = np.zeros(nDay, dtype=int)  # initialize the counts per fraction of a day
    dT = rs.emjd - mjdRef   # get fraction of day of this event
    dT = dT * nDay          # get index to counts of that time
    iT = int ( dT)          # convert to integer
    countsInDt[iT] = 1      # flag that event
    aDir = { 'dir': "pi1-events-25Dec13",    # name of directory
             'events' : anEvent,   # names of events in directories
             'n': 1,               # count of events in directory
             'mjds': mjds,            # array of MJDs of these events
             'counts' : countsInDt
            }
    # pass flash count and directory name of events
    eventDirs = { 0 : aDir }         # create 3 identical events
    mjds = [ rs.emjd, mjdRef + .1]
    iT = int( .1 * nDay)
    countsInDtB = copy.deepcopy(countsInDt)
    countsInDtB[iT] = 1
    bDir  = { 'dir': "pi2-events-25Dec13",    # name of directory
             'events' : anEvent,   # names of events in directories
             'n': 2,               # count of events in directory
             'mjds': mjds,            # array of MJDs of these events
             'counts' : countsInDtB
            }
    eventDirs.update({ 1 : bDir })
    mjds = [ rs.emjd, mjdRef + .1, mjdRef + .2]
    iT = int( .2 * nDay)
    countsInDtC = copy.deepcopy(countsInDtB)
    countsInDtC[iT] = 1
    cDir  = { 'dir': "pi3-events-25Dec13",    # name of directory
             'events' : anEvent,   # names of events in directories
             'n': 3,               # count of events in directory
             'mjds': mjds,            # array of MJDs of these events
             'counts' : countsInDtC
            }
    eventDirs.update({ 2 : cDir })
    
    logFile, directoryDate = logDirectory( calendar, outDirName, eventLogName)

    raTransit = ((42./60) + (44./3600.))*15
    decTransit = 41.  + (16./60.) + (8./3600)
    transitObject ="M31"
    plotHistogram( directoryDate, calendar, nDir, rs, nDay, \
                   mjdRef, transitObject, raTransit, decTransit, eventDirs,
                   1, [ rs.emjd], [ 1], [ rs.gallon ], [rs.gallat], \
                   [" ", " ", " "], 1)
    logFile.close()

# now count all events happening within .1 degrees of other events.

if __name__ == "__main__":
    main()

