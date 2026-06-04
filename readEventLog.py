#/bin/python
# read in log files from matchevents.py and plot and analyze
#HISTORY
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

# locally defined class and functions
import GalacticHistogram
import radioastronomy
import nameToIndex

# temporarily turn on/off debugging
verbose = True
verbose = False

nnplot = 0
MAXDIR = 20

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

def logSummary( logFileName):
    """
    Check if the log File is present and return a set of arrays of with length
    equal to the number of days in the log file
    """

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

    # count lines with dates
    nDay = 0
    dates = []
    teles = []
    events = []
    tens   = []
    hundreds = []
    thousands = []
    lastDate = ""
    for aline in lines:
        # get rid of extra spaces at beginning, end and middle of string
        oneline = " ".join(aline.split())
        if oneline[0] == '#' and oneline[1] == 'F':
            continue
        elif oneline[0] == "#":
            parts = oneline.split(" ")
            nparts = len(parts)
            if nparts != 6:
                print("Do not recognize this line: %s (%d)" % (oneline, nparts))
                continue
            # get #26May01 part
            adate = parts[0]
            adate = adate[1:]  # skip #
            nDay = nDay + 1
            dates.append( adate)
            teles.append( int(parts[1]))
            events.append( int(parts[2]))
            tens.append( int(parts[3]))
            hundreds.append( int(parts[4]))
            thousands.append( int(parts[5]))
            if adate == lastDate:
                print("Warning Date is duplicated in Event List!: %s" % (date))
            lastDate = adate
            
    # end of all lines
    # now convert lists to arrays
    teles = np.array( teles)
    events = np.array( events)
    tens = np.array( tens)
    hundreds = np.array( hundreds)
    thousands = np.array( thousands)

    if nDay < 1:
        print("No Events in log File, Exiting!")
        exit()
    else:
        print("Read %d Days in Log File" % (nDay))

    return dates, teles, events, tens, hundreds, thousands
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

def plotHistogram( nDay, dates, teles, events, tens, hundreds, thousands):
    """
    plot several histograms of the event count versus time of day
    where:
    nDay = number of days of event summary
    dates - YYmmmDD string for each date of observation
    
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    global calendar
    fig, ax = plt.subplots( figsize=(12,6))

#    barcolors = ["lightcyan", "wheat", "greenyellow", "mistyrose", \
#                 "wheat", "lightgrey", "mintcream"]
    barcolors = ["gold", "lime", "dodgerblue", "darkmagenta", \
                 "mistyrose", "wheat", "lightgrey", "mintcream"]

    xaxis = np.arange(nDay) + 1
    plt.bar( xaxis, events, label="Events", color=barcolors[0])
    plt.bar( xaxis, tens, label="N > 10", color=barcolors[1])
    plt.bar( xaxis, hundreds, label="N > 100", color=barcolors[2])
    plt.bar( xaxis, thousands, label="N > 1000", color=barcolors[3])
    plt.bar( xaxis, teles, label="Telescopes", color='none', edgecolor='black')

    plt.legend()
    plt.title("%s to %s   Event Counts" % (dates[0], dates[nDay-1]),fontsize=18)
    plt.xlabel("Day", fontsize=18)
    plt.ylabel("Number of Events", fontsize=18)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.xlim(.5, nDay+.5)
#    xticks = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
#    plt.xticks( float( np.array(xticks)), str(np.array(xticks)))
#    nticks = len(xticks)
#    strticks = xticks
#    for i in range(nticks):
#        strticks[i] = str( xticks[i])
#        xticks[i] = float(xticks[i])
#    ax.set_xticks( xticks)

    #    ax.set_xticklabels( strticks)
#    ax.tick_params(which='major', width=3)
#    ax.tick_params(which='minor', width=1)
#    ax.tick_params(which='major', length=5)
#    ax.tick_params(which='minor', length=3, color='b')
#    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.25))
#    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
#    fig = plt.gcf()
#    fig.savefig(plotfile, bbox_inches='tight')
#    print("Saved match summary plot: %s" % (plotfile))
    plt.show()
    return
# end of plotHistogram

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

        parser.add_argument('-v', '--verbose', default=False,
                            action='store_true')

        # end of readArgParse()
        return parser.parse_args()

    inputs = readArgParse()

    nLogs = 0
    if inputs.Logs != None:
        eventLog = inputs.Logs
        nLogs = len(eventLog)

    if nLogs < 1:
        print("No event log names supplied!")
        exit()

    for iLog in range(nLogs):
        aLog = eventLog[iLog]
        dates, teles, events, tens, hundreds, thousands = logSummary(aLog)
        nDay = len( dates)
        plotHistogram( nDay, dates, teles, events, tens, hundreds, thousands)
        
    # end of main()
    
# if testing reading of event logs
if __name__ == "__main__":
    main()

