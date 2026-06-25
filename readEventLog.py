#/bin/python
# read in log files from matchevents.py and plot and analyze
#HISTORY
#26Jun16 GIL Event counts sorted, now overplot
#26Jun15 GIL sort all types of events F=<10, X>=10, C>=100, M>=1000
#26Jun08 GIL clean up event type matching, add labels
#26Jun07 GIL select particular events for histogram
#26Jun04 GIL Plot events as a function of galactic latitude
#26Jun04 GIL plot events per day
#26Jun03 GIL read only the dates summaries
#26Jun02 GIL Initial version from matchEvent
#

import sys
import os
import datetime
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.coordinates import AltAz
import argparse
import radioastronomy
from jdutil import date_to_mjd
from findDate import findDirs, readEventsInDir
import utcOffset

# locally deftned class and functions
import GalacticHistogram
import radioastronomy
import nameToIndex

# temporarily turn on/off debugging
verbose = True
verbose = False
barcolors = ["gold", "orange", "magenta", "crimson", "darkmagenta", \
                 "mistyrose", "wheat", "lightgrey", "mintcream"]
barlabels = ["Events", "N < 10", "N > 9", "N > 99", "N > 999", \
                 "mistyrose", "wheat", "lightgrey", "mintcream"]
nnplot = 0
MAXDIR = 20
#26Jun15 GIL sort all types of events 0=<9, X>=10, C>=100, M>=1000
0  # Indicies for a few events
X=1  # 10s of events
C=2  # 100s of events
M=3  # 1000s of events
F=4  # Fit to events
#There are 3 types of lines in an eventLog.
#1 The Date Summary, with 5 entries
#
"""
#26May01  9     3     0     0     0
26May01 pi11-events-26May01/26-05-01T033831_842.eve   61161.1518  3 ' '  121.37   78.42  193.26   38.71
26May01 pi11-events-26May01/26-05-01T224805_327.eve   61161.9500  4 '2'  182.11   30.21  121.19   38.57
26May01 pi9-events-26May01/26-05-01T235119_131.eve   61161.9949  4 '5'  184.11   42.51  137.08   38.64
#26May02  9    15     4     5     0
#FIT 26May02  12.357+/-0.245  0.576 122125.015       61162.5149  1526.7+/-154.5
26May02 pi11-events-26May02/26-05-02T001854_657.eve   61162.0132  3 '2'  184.25   47.92  144.02   38.66
26May02 pi11-events-26May02/26-05-02T010740_036.eve   61162.0470  5 '2'  182.75   57.45  156.28   38.69
26May02 pi7-events-26May02/26-05-02T102930_279.eve   61162.4372  3 ' '   73.12    6.13  297.43   38.29
"""
def parseDateSummary( oneline, lastDate):
    """
    parse line with format like:
    """
    parts = oneline.split(" ")
    nparts = len(parts)
    if nparts != 6:
        print("Do not recognize this line: %s (%d)" % (oneline, nparts))
        adate = ""
        # bad line, no counts
        counts = [ 0, 0, 0, 0, 0, 0]
    else:
        # get #26May01 part
        adate = parts[0]
        adate = adate[1:]  # skip #
        teles = int(parts[1])
        events = int(parts[2])
        tens = int(parts[3])
        hundreds = int(parts[4])
        thousands = int(parts[5])
        # single event count is total minus all others
        ones = events - (tens + hundreds + thousands)
        if adate == lastDate:
            print("Warning Date is duplicated in Event List!: %s" % (date))
        lastDate = adate
        counts = [teles, events, ones, tens, hundreds, thousands]
    # end of parseDateSummary
    return adate, lastDate, counts

def logSummary( logFileName, flag = " "):
    """
    Check if the log File is present and return a set of arrays of with length
    equal to the number of days in the log file
    Optionally select event types:
    flag = 'F'  - Fit events
    flag = ' '  - all events
    flag = '1'  - all integer flags
    flag = 'X'  - all events with 10 to 99 flashes
    flag = 'C'  - all events with 100 to 999 flashes
    flag = 'M'  - all events with 1000 or more flashes
    """

    if flag == ' ':
        print("Selecting all Events")
    else:
        print("Selecting only events of type '%s'" % (flag))
    
    logOK = False   # assume can not read log file
    if os.path.exists(logFileName):
        if os.path.isfile(logFileName):
            try:
                logFile = open(logFileName, 'r')
                logOK = True
            except:
                logOK = False
    if logOK:
        print("Open Log file: %s" % (logFileName))
    else:
        print("Can not Read Log file: %s" % (logFileName))
        return -1
    # now read all the lines in the file(s)
    lines = logFile.readlines()
    logFile.close()

    # count lines with dates  Type A lines
    nDay = 0
    dates = []
    counts = [ 0, 0, 0, 0, 0]
    eventCounts = []
    lastDate = ""
    # type B lines = individual events
    bdates = []
    bfiles = []
    bmjds  = []
    bteles = []
    bflags = []
    blons  = []
    blats  = []
    bras   = []
    bdecs  = []
    nB = 0
    # type C lines == FITS
    cdates = []
    ctimes = []
    ctwidths = []
    cdurs  = []
    ctods  = []
    cmjds  = []
    cnevents = []
    cnrms  = []
    nC = 0
    nPrint = 0
    nFound = 0
    nSkip  = 0
    # now identify and parse each line
    for aline in lines:
        # get rid of extra spaces at beginning, end and middle of string
        inline = aline
        oneline = " ".join(inline.split())
        # if a #FIT line        
        if oneline[0] == '#' and oneline[1] == 'F':
            nC = nC + 1
            continue
        # if summary of events on a date
        elif oneline[0] == "#":
            adate, lastDate, counts = parseDateSummary( oneline, lastDate)
            if adate != "":
                dates.append( adate)
                eventCounts.append( counts)
                nDay = nDay + 1
        elif oneline[0] == "2":  # If line with a date line. ie: 26Jun04
            parts = oneline.split(" ")
            nparts = len(parts)
            if nparts != 9:
                print("Line not recognized: ")
                print("%s" % (oneline))
                print("%s" % (aline))
                continue
            else:
                aFlag = parts[4]   # expect flag to have format 'X'
                aFlag = str(aFlag[1])   # get only single character, X
                #print( "Comparing flag %s to %s" % (aFlag, flag))
                # seeking one flag of a specific type and no match
                matchFlag = False
                if flag == " ":
                    matchFlag = True
                elif flag == "1":
                    matchFlag = aFlag in "123456789"
                elif flag == aFlag:
                    matchFlag = True
                
                nPrint = nPrint + 1
                if matchFlag:
                    if nFound < 1:
                        print( "   Found flag '%s' == '%s'" % (aFlag, flag))
                        nFound = 1
                else:
                    if nSkip < 1:
                        print( "Skipping flag '%s' != '%s'" % (aFlag, flag))
                        nSkip = 1
                if matchFlag == False:
                    continue

                bdates.append( parts[0])
                bfiles.append( parts[1])
                bmjds.append( float( parts[2]))
                bteles.append( int(parts[3]))
                bflags.append( parts[4])
                blons.append( float( parts[5]))
                blats.append( float( parts[6]))
                bras.append( float( parts[7]))
                bdecs.append( float( parts[8]))
                nB = nB + 1
                
    # end of all lines
    
    # now conver B lists to arrays
    bmjds  = np.array( bmjds)
    bteles = np.array( bteles)
    blons  = np.array( blons)    
    blats  = np.array( blats)    
    bras   = np.array( bras)    
    bdecs  = np.array( bdecs)    
    if nDay < 1:
        print("No Days of Events in log File, Exiting!")
        exit()
    else:
        if flag == " ":
            print("Read %d Days in Log File" % (nDay))

    if nB > 0:
        print( "Found %d events of type '%s'" % (nB, flag))
    else:
        print( "Found No events of type '%s'" % (flag))
        exit()
        
    return dates, eventCounts, \
        bdates, bfiles, bmjds, bteles, bflags, blons, \
        blats, bras, bdecs

    # end of logSummary()
    
def logEvent( logFile, files, event, gallon, gallat):
    """
    logEvent() adds one event to the event log
    where
    logFile = file pointer to open log file
    files   = composed of this event
    event   = event to be added to the log
    """

    """
        print("#FIT %s %7.3f+/-%.3f %6.3f %s     %12.4f %7.1f+/-%.1f" % \
              (utcYYparts[0], time*24., stdDev*24.,\
               fwhm*24., utcYYparts[1], mjd, peak, rms), file=logFile)
    """
    """
    print( "#%s %2d %5d %5d %5d %5d" % \
           (calendar, nDir, nGroup, nTen, nHundred, nThousand), \
           file=logFile)
    """
    """
    print("%s%s" % (aspace, afile), file=logFile)
    """
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

def plotHistogram( nDay, dates, eventCounts, title = ""):
    """
    plot several histograms of the event count versus day
    where:
    nDay = number of days of event summary
    dates - YYmmmDD string for each date of observation
    eventCounts = lists of 6 integers: telescopes, events, ones, tens,
       hundreds, thousands)
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    fig, ax = plt.subplots( figsize=(8,6))

    teles = [0] * nDay            # make arrays to fill with counts
    events = [0] * nDay
    ones = [0] * nDay
    tens = [0] * nDay
    hundreds = [0] * nDay
    thousands = [0] * nDay
    #### 
    xaxis = np.arange(nDay) + 1
    for i in range( nDay):
        teles[i] = eventCounts[i][0]
        events[i] = eventCounts[i][1]
        ones[i] = eventCounts[i][2]
        tens[i] = eventCounts[i][3]
        hundreds[i] = eventCounts[i][4]
        thousands[i] = eventCounts[i][5]
    plt.bar( xaxis, events, label=barlabels[0], color=barcolors[0])
    plt.bar( xaxis, ones, label=barlabels[1], color=barcolors[1])
    plt.bar( xaxis, tens, label=barlabels[2], color=barcolors[2])
    plt.bar( xaxis, hundreds, label=barlabels[3], color=barcolors[3])
    plt.bar( xaxis, thousands, label=barlabels[4], color=barcolors[4])
    plt.bar( xaxis, teles, label="Telescopes", color='none', edgecolor='black')

    plt.legend()
    if title != "":
        plt.title(title, fountsize=18)
    else:
        plt.title("%s to %s Counts" % \
                  (dates[0], dates[nDay-1]),fontsize=18)
    plt.xlabel("Day", fontsize=18)
    plt.ylabel("Number of Events", fontsize=18)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.xlim(.5, nDay+.5)
    plt.show()
    return
# end of plotHistogram

def plotGalactic( bdates, bfiles, bmjds, bteles, bflags, blons, \
                  blats, bras, bdecs, flag = " ", title = "", date="", \
                  color="", label="", degrees = 4., show=True):
    """
    plot Events in a histograms of Galactic latitude bin
    where:
    nDay = number of days of event summary
    dates - YYmmmDD string for each date of observation
    
    """
    glh = GalacticHistogram.GalacticLatitudeHistogram( bin_width_deg=degrees)

    nEvents = len( bdates)
    for iEvent in range( nEvents):
        glh.add_measurement( blats[iEvent], 1.)

    print("Ploting events of type: %s" % (label))
    
    centers, allEvents = glh.plot( plotType = GalacticHistogram.EVETYPE, \
                                   date=date, color=color, \
                                   title=title, label=label, show=show)

    allEvents = np.array( allEvents)
    centers = np.array( centers)
    return centers, allEvents
# end of plotGalactic

def plotEventHistogram( latitudes, scaleEvents, flags, width=4., title = ""):
    """
    plot all Events Types on a single plot, versus galactic latitude
    where:
    latitude = latitude bins
    scaleEvents = two dimensional set of of number of events in each range
    flags = character specifying event type
    title = additional title string, usually with dates
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    fig, ax = plt.subplots( figsize=(8,6))

    events = scaleEvents[0]
    ones = scaleEvents[1]
    tens = scaleEvents[2]
    hundreds = scaleEvents[3]
    thousands = scaleEvents[4]
    # If your data has an offset/mean

    ####
    xaxis = latitudes
    plt.bar( xaxis, events, label="Events", width=width, color=barcolors[0])
    plt.bar( xaxis, ones, label="N < 10", width=width, color=barcolors[1])
    plt.bar( xaxis, tens, label="N > 9", width=width, color=barcolors[2])
    plt.bar( xaxis, hundreds, label="N > 99", width=width, color=barcolors[3])
    plt.bar( xaxis, thousands, label="N > 999", width=width,color=barcolors[4])

    n = len(xaxis)
    xmin = 999
    xmax = -999
    imin = 0
    imax = 0
    # now find ranges of longitude where data are present
    for i in range(n):
        if events[i] > 0:
            if xaxis[i] < xmin:
                xmin = xaxis[i]
                imin = i
            if xaxis[i] > xmax:
                xmax = xaxis[i]
                imax = i
    plt.xlim( xmin-(width/2.), xmax+(width/2.))
    plt.legend()
    fulltitle = "Event Rate vs Galactic Latitude"
    if title != "":
        fulltitle = fulltitle + " : " + title
    plt.title(fulltitle,fontsize=18)
    plt.xlabel("Galactic Latitude (deg)", fontsize=18)
    plt.ylabel("Scaled Number of Events", fontsize=18)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.show()
    # end of plotEventHistogram()
    for i in range( 5):
        data = np.array(scaleEvents[i])
        data = data[imin:(imax+1)]
        xs  =  xaxis[imin:(imax+1)]
        mean = np.mean( data)
        rms = np.sqrt(np.mean(data**2))

        # Mathematically: Variance = RMS^2 - Mean^2
        sigma = np.sqrt(rms**2 - mean**2)
        print("%8s: Average Counts: %6.2f +/- %5.2f" % \
              (barlabels[i], mean, sigma))
        idatamax = np.argmax(data)
        datamax = data[idatamax]
        latmax = xs[idatamax]
        nsigma = (datamax - mean)/sigma
        print("%8s  Maximum offset: %6.2f (%.2f simga) at %6.2f deg" % \
              (" ", datamax, nsigma, latmax))
        # end for all flags, compute sigma
    return

def main():
    '''
    main() function plot summaries of event logs
    '''
    import argparse
    import textwrap

    epilogText = textwrap.dedent('''Read all matched events in log files and plot histograms of events versus Galactic Latitude and other parameters

    Glen Langston 2026 June 2

Example:
    python readEventLog.py -v /media/karl/daily/2026/26May/26May-Events.txt
    ''')

    def readArgParse():
        parser = argparse.ArgumentParser(
            description="Read set of event summary logs and plots histograms of data versus different parameters.  The main test is versus Galactic Latitude to test whether the events are likely of Galactic Origin.",
            epilog=epilogText,
            formatter_class=argparse.RawDescriptionHelpFormatter
            )

        # Example of another optional string
        parser.add_argument(
            "--title",
            type=str,
            help="Optional title for top of output plots"
        )

        parser.add_argument("--hist", type=str,
                            help="Name of file containing the histogram of event durations during the times of event detections.   Used to normalize the Galactic Latitude dependence of the event rates."
                            )

        parser.add_argument('Logs', nargs='*', default=[""], 
                            type=str,
                            help="List of log files to evaluate"
                            )

        parser.add_argument('-f', '--flag', default=" ", 
                            type=str,
                            help="Select type of event Flag to plot.  Options are 1, X, C, M.   Default is all events"
                            )

        parser.add_argument('-v', '--verbose', default=False,
                            action='store_true')

        # end of readArgParse()
        return parser.parse_args()

    inputs = readArgParse()

    # select type of events to plot
    flag = inputs.flag
    if flag == " ":
        print( "Selecting all events for histogram")
    else:
        print( "Selecting events with flag: %s" % (flag))

    nLogs = 0
    if inputs.Logs != None:
        eventLog = inputs.Logs
        nLogs = len(eventLog)

    if nLogs < 1:
        print("No event log names supplied!")
        exit()

    # if correcting counts for time in each latitude range
    if inputs.hist != None:
        fullPath = inputs.hist
        print("Duration Histogram: %s" % (fullPath))
        # read time in each galactic latitude range
        date, dataDir, telIndex, degree, centers, samples = \
            GalacticHistogram.readHistogram( fullPath)
    nCenters = len( centers)

    # create a set of arrays for each type
    eventTypes = []
    # for each type of event, " " == all events
    flagList = [" ", "1", "X", "C", "M"]
    iFlag = 0
    for aFlag in flagList:
        # restart the gathering of events to plot
        adates = []
        files = []
        mjds = []
        teles = []
        flags = []
        lons = []
        lats = []
        ras = []
        decs = []
        
        # read all log files 
        for iLog in range(nLogs):
            aLog = eventLog[iLog]
            dates, eventCounts, \
                bdates, bfiles, bmjds, bteles, bflags, blons, \
                blats, bras, bdecs = logSummary(aLog, aFlag)
            nDay = len( dates)
            if aFlag == " ": 
                plotHistogram( nDay, dates, eventCounts)

            nType = len(bdates)
            for iType in range( nType):
                adates.append( bdates[iType])
                files.append( bfiles[iType])
                mjds.append( bmjds[iType])
                teles.append( bteles[iType])
                flags.append( bflags[iType])
                lons.append( blons[iType])
                lats.append( blats[iType])
                ras.append( bras[iType])
                decs.append( bdecs[iType])
        # end for all types of events 
        # now all events of this type are gathered, plot them
        title = "Galactic Latitude Histogram: %s" % (barlabels[iFlag])
        date = "%s - %s" % (dates[0], dates[nDay-1])
        
        latitudes, allEvents = plotGalactic( adates, files, \
                                             mjds, teles, bflags, \
                                             lons, lats, ras, decs, \
                                             flag=aFlag, date=date, \
                                             title=title, \
                                             color=barcolors[iFlag], \
                                             label=barlabels[iFlag],
                                             show = False)

        eventTypes.append( allEvents)
        iFlag = iFlag + 1
        # end all types of events
    nLatitudes = len( latitudes)
    if (nCenters != nLatitudes):
        print( "Miss-match: spectra and event latitudes: %d != %d" \
               % (nCenters, nLatitudes))
        exit()
    maxSamples = max( samples)

    scaleEvents = []
    iFlag = 0
    # now for each type of flag, scale by time in range of galactic latitude
    for aFlag in flagList:
        # for all latitude bins
        calEvents = np.zeros( nCenters)
        for i in range(nCenters):
            allEvents = eventTypes[iFlag]
            if samples[i] > 0.:
                calEvents[i] = maxSamples * allEvents[i] / samples[i]
        scaleEvents.append( calEvents)
        iFlag = iFlag + 1
    # now plot normalized event rates per latitude

    title = "%s - %s" % (dates[0], dates[nDay-1])
    plotEventHistogram( latitudes, scaleEvents, flags, title=title)
    # end of main()
    
# if testing reading of event logs
if __name__ == "__main__":
    main()

