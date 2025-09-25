#Python find matchs in data directories
#HISTORY
#25Aug11 GIL test for single events, to avoid len(1) failing
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

MAXDIR = 20

def findDirs( dirs, calendar):
    """
    findDirs() finds the names of directories matching the calendar flag
    """
    dirs = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""]
    nDir = 0
    if calendar == "":
        try:
            nDir = len(dirs)
        except:
            dirs[0] = dirs
            nDir = 1
        return dirs, nDir
    else:
        dirs = []
    # if a date name, createin directories
        for idir in range(1, 13):
            adir = "pi%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs.append(adir)
                nDir = nDir + 1
                print("Found directory: %s" % adir)
            if nDir >= MAXDIR:
                break

        for idir in range(1, 15):
            adir = "odroid%d-events-" % idir
            adir = adir + calendar
            if os.path.isdir(adir):
                dirs.append(adir)
                nDir = nDir + 1
                print("Found directory: %s" % adir)
            if nDir >= MAXDIR:
                break

    if nDir < 2:
        print("Can not match events in less than 2 directories")
# end of findDirs
    return dirs, nDir

def fileNameToMjd( filename):
    """
    Extract the event time from the file name
    """
    
# expecting a name like pi6-events-25Aug08/25-08-09T013527_390.eve
    fileparts = filename.split("/")
    nparts = len(fileparts)
    # now have only 25-08-09T013527_390.eve
    datetime = fileparts[ nparts-1]
    # split date and time
    dateparts = datetime.split("T")
    date = dateparts[0]
    # now have 25-08-09
    ymdparts = date.split("-")
    yyy = int(ymdparts[0])
    mmm = int(ymdparts[1])
    ddd = float(ymdparts[2])
#    print( ymdparts, yyy, mmm, ddd)
    # next work on 013527_390.eve
    timeparts = dateparts[1].split(".")
    time = timeparts[0]
    # time is now 013527_390
    secondparts = time.split("_")
    if not len(secondparts) == 2:
        print( "time not what is expected %s ", (dateparts[1]))
        millisec = 0.
    else:
        millisec = float(secondparts[1])
    # now work on hhmmss
    hhmmss = secondparts[0]
    hh = float(hhmmss[0:2])
    mm = float(hhmmss[2:4])
    ss = float(hhmmss[4:6])
    seconds = (hh*3600.) + (mm*60.) + ss + (millisec/1000.)
    dayfraction = seconds / 86400.
    mjd = date_to_mjd( 2000.+yyy, mmm, ddd + dayfraction) 
#    print( filename, hh, mm, ss, millisec, mjd)
# end of filename to mjd
    return( mjd)

def readEventsInDir( directory, verbose = False):
    """
   Read in all events in a file list and return arrays
    """
#   only read the az, el once for each directory
    readAzEl = True
    
#   only examine event files
    events = list(glob.glob(os.path.join(directory,'*.eve')))
    try:
        nEve = len(events)
    except:
        nEve = 0
        
    if nEve < 1:
        print( "No Events in Directory; %s" % (directory))
        return nEve, "", 0.
    else:
        if verbose:
            print("Found %5d events in directory: %s" % (nEve, directory))

    mjds = np.zeros(nEve)

    kkk = 0
    for filename in events:
        fullname = filename
        if verbose and kkk < 3:
            print("%5d: %s" % (kkk, fullname))
        if readAzEl:
            rs = radioastronomy.Spectrum()
            rs.read_spec_ast(fullname)
            readAzEl = False
# to speed up processing, only use file name to get date of event.
        mjd = fileNameToMjd( fullname)
        mjds[kkk] = mjd
        if (100 * int(kkk/100) == kkk) and verbose:
            print("Event %5d: %12.9f: %7.3f+/-%5.3f" % \
                  (kkk, rs.emjd, rs.epeak, rs.erms))
        kkk = kkk + 1
    mjds = mjds[0:kkk]
    nEve = kkk
    if verbose:
        if nEve < 1:
            print("No Events found")
        else:
            print("N events: %d, last file %s, last mjd: %f" % ( nEve, events[nEve-1], mjds[nEve-1]))
    return nEve, events, mjds
# end of readEventsInDir

