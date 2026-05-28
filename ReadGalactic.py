#/bin/python
#Calculate histogram of galactic latitude ranges, weighted by the time in range
#optionally plot histogram
#to test use:  python ReadGalactic.py data
#HISTORY
#26May28 Gil record peak intensity for duration of events
#26May27 Gil finish adding argparse, compare and write histograms
#26May26 Gil add date annotate the plot with telescope and elevation
#26May25 Gil annotate the plot with telescope and elevation
#26May24 Gil initial version
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import argparse
import astropy

# locally defined class and functions
import GalacticHistogram
import radioastronomy
import nameToIndex

def getDirNameParts( dirname):
    """
    Assuming a standard name format, including telescope and date
    return these parts of the directory name
    """
    telIndex = nameToIndex.getIndex( dirname)

    # get month of observation for directory with format like:
    # /media/karl/pi11-data-26Mar26
    parts = dirname.split("-")
    if len(parts) < 2:
        date = ""
        print("No date in directory name")
    else:
        # else date is last segment of name
        date = parts[len(parts)-1]

    return date, telIndex

#Read directory of observations and compute histogram of spectra/events
def ReadGalactic( dirname, degrees=4., show=True):
#Optionally run a diagnotic test
    """
    Main executable testing histogram
    """
    
    date, telIndex = getDirNameParts( dirname)
        
    # creast a list of files in the input directory
    from os import listdir
    from os.path import isfile, join
    files1 = [f for f in listdir(dirname) if isfile(join(dirname, f))]
    
    print("Degree width of histogram bins: %0.1f " % (degrees))

    print("%5d Files in Directory: %s" % (len(files1), dirname))

    # initialize function to read spectra, events and data header
    rs = radioastronomy.Spectrum()

    obs1s = files1
    nFiles = len( files1)
    types = np.zeros( nFiles+1)
    nFile = 0
    nHot = 0
    minel = 99.
    maxel = -99.
    minaz = 365.
    maxaz = -365.
    # now for all files
    for filename in files1:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        isObs = isEve = isHot = False
        # count observations and
        if parts[nparts-1] == "ast":
            isObs = True
            types[nFile] = GalacticHistogram.ASTTYPE
        elif parts[nparts-1] == "eve":
            isEve = True
            types[nFile] = GalacticHistogram.EVETYPE
        elif parts[nparts-1] == "hot":
            isHot = True
            nHot = nHot + 1
            types[nFile] = GalacticHistogram.HOTTYPE
        else:
            continue
        # if here, then the file is one of known types
        obs1s[nFile] = filename
        # count only astronomical Types found
        nFile = nFile + 1

    # start a new histogram for this directory
    glh = GalacticHistogram.GalacticLatitudeHistogram( bin_width_deg=degrees)
    # either count spectral observations or now get event signficances
    iii = 0
    total = 0.    # time observing spectra
    tcal  = 0.    # time calibrating spectra
    etotal = 0.    # time observing spectra
    # for all fies in this directory
    for filename in obs1s[0:nFile]:
        fullname = join(dirname, filename)
        rs.read_spec_ast( fullname, headerOnly=True)
        rs.azel2radec()    # compute ra,dec from az,el
        if 50 * int(iii/50) == iii:
            print("File %5d: az,el = %5.1f,%6.1f => galactic %5.1f,%6.1f" % \
                  (iii, rs.telaz, rs.telel, rs.gallon, rs.gallat))
        tint = rs.durationSec
        lat = rs.gallat
        if types[iii] == GalacticHistogram.ASTTYPE:
            total = total + tint
            # now add to histogram 
            glh.add_measurement( lat, tint)
        elif types[iii] == GalacticHistogram.HOTTYPE:
            tcal = tcal + tint
        elif types[iii] == GalacticHistogram.EVETYPE:
            etotal = etotal + rs.epeak
            if rs.erms > 0.:
                ePeak = rs.epeak/rs.erms    # durations are nearly zero, use SNR
            else:
                print("Events: %7.3f +/- %.3f" % (rs.epeak, rs.erms))
                ePeak = rs.epeak    # computation error, use peak intensity
            # integrate SNRs
            etotal = etotal + ePeak
            # now add to histogram 
            glh.add_measurement( lat, ePeak)
        iii = iii + 1  # step through files to examime

        # only keep track of min and max elevation if not calibrating
        if types[iii] == GalacticHistogram.HOTTYPE:
            continue
        # keep track of location of observations
        if rs.telel > maxel:
            maxel = rs.telel
        if rs.telel < minel:
            minel = rs.telel
        if rs.telaz > maxaz:
            maxaz = rs.telaz
        if rs.telaz < minaz:
            minaz = rs.telaz

    if nFile < 1:
        print(""")
        print("Not enough data for a histogram, exiting...")
        print(""")
        exit()

    # must assume all inputs are of the same type
    plotType = types[0]
    if plotType == GalacticHistogram.ASTTYPE:
        title = "Tel %2d: Galactic Latitude Histogram" % telIndex
    elif plotType == GalacticHistogram.EVETYPE:
        title = "Tel %2d: Galactic Event Histogram" % telIndex
    else:
        title = ""
    if maxaz == minaz:
        title = title + (" Az: %.1f" % maxaz)
    else:
        title = title + (" Az: %.1f to %.1f" % (minaz, maxaz))
    if maxel == minel:
        title = title + (" El: %.1f" % maxel)
    else:
        title = title + (" El: %.1f to %.1f" % (minel, maxel))
    # now plot histogram
    centers, norm_hist = glh.plot( date=date, title=title, plotType=plotType)

    print(title)
    print("Fraction of time Observing %.1f %s" % \
          (glh.total_time_sec/864., "%"))
    print("Fraction of Calibration    %.1f %s" % \
          (tcal/864., "%"))
    if plotType == GalacticHistogram.ASTTYPE:
        norm_hist = glh.normalized() * 1440.
    else:
        norm_hist = glh.normalized()
    return date, telIndex, centers, norm_hist

if __name__ == "__main__":
    '''
    main() function defined to test counting 
    '''
    import argparse
    import textwrap

    epilogText = textwrap.dedent('''Compute the Galactic Latitude Distribution of Events detected during NSF Watch observations.

    Two sets of observations are compared.
    1) The Spectral Line observations, which should always be on-going, so are a good determination that the telescope is functioning.
    2) The Events detected by this telescope during the same date.

    Glen Langston 2026 May 26

Example:
    python ReadGalactic.py -d 4. /media/karl/pi1-data-26May25 /media/karl/pi1-events-26May25
    ''')

    def readArgParse():
        parser = argparse.ArgumentParser(
            description="Read a directory of Spectra and a Directory of Events and deduce the Galactic Latitude distribution of the events, including a calibration of the relative time observing different Galactic Latitudes.",
            epilog=epilogText,
            formatter_class=argparse.RawDescriptionHelpFormatter
            )

        # Optional distance to source, for calculation of molecule mass
        parser.add_argument(
            "-l", "--log",
            type=str,
            default = "/home/karl/daily",
            help="Output directory root for histograms"
        )        

        # Optional distance to source, for calculation of molecule mass
        parser.add_argument(
            "-d", "--degrees",
            type=str,
            default = "4.", 
            help="Optionally provide galactic Latitude Bin Width (degrees)"
        )        

        # Example of another optional string
        parser.add_argument(
            "--title",
            type=str,
            help="Optional title for top of output plots"
        )

        parser.add_argument('SpectralDir', type=str,
                            help="Name of directory containing spectra for the time of event detection.  Used to calibrate Galactic Latitude distribution of Events"
                            )

        parser.add_argument('EventDirs', nargs='*', default=[""], 
                            type=str,
                            help="List of direcetories containing event detections"
                            )

        parser.add_argument('-v', '--verbose', default=False,
                            action='store_true')
        
        return parser.parse_args()

    args = readArgParse()

    degrees = float( args.degrees)

    # define
    logDirectory = args.log
    
    # first read spectra to normalize event counts vs Galactic Latitude
    if args.SpectralDir != None:
        spectra = str( args.SpectralDir)
        # from spectra data name deduce date and telescope index
        date, telIndex = getDirNameParts( spectra)
        # create full path to possibly existing histogram log file
        fullPath = GalacticHistogram.createLogName( logDirectory, date, telIndex)
        filePath = Path(fullPath)
        if filePath.exists():
            print("Reading Previous Spectra Histogram")
            date, dataDir, telIndex, degree, centers, samples = \
                GalacticHistogram.readHistogram( fullPath)
        else:
            date, telIndex, centers, samples = \
                ReadGalactic(spectra, degrees=degrees)

            # if data read OK, then write histogram
            if len(date) < 2:
                print("Date not found in name, not logging: %s " % (date))
            else:
                GalacticHistogram.writeHistogram( logDirectory, \
                                                  date, spectra, \
                                                  telIndex, degrees, \
                                                  centers, samples)

    nEvent = 0
    if args.EventDirs != None:
        eventDirs = args.EventDirs
        for anEventDir in eventDirs:
            eventDate, eventIndex, eventCenters, eventSamples = \
                ReadGalactic(anEventDir, degrees=degrees)
            if nEvent == 0:
                eventSums = eventSamples
            else:
                eventSums = eventSums + eventSamples
            nEvent = nEvent + 1

            

