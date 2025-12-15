#Python find matchs in data directories
#HISTORY
#25Sep24 GIL Fix Not verbose only functions
#25Sep23 GIL initial version from code in matchevents.py

import sys
import shutil
import os
import copy
import numpy as np
from findTime import *

# time offsets are in units of days.
OneMjdSec = float(1./86400.)
NOMATCH = -9999
MAXDT = 365.   # 1 year

def findpairs( EventDir1, EventDir2, tOffset = OneMjdSec):
    """
    findpairs() - match events in pairs of event structures
    Inputs are Python data structures describing the events in the directories
    Outputs are indicies of matches to 2nd directory and time offsets between
    """
    # extract number of events in each directory
    n1 = EventDir1['n']
    n2 = EventDir2['n']
    # make sure there are plenty of indicies, set to no-match
    ii12s = np.zeros(n1, dtype=int) + NOMATCH
    # init minimum time offsets   ; One index for every entry in directory 1
    dt12s = np.zeros(n1) + MAXDT

    mjd1s = EventDir1['mjds']
    mjd2s = EventDir2['mjds']

    if n1 > 1:
        isArray1 = True
    else:
        isArray1 = False
    if n2 > 1:
        isArray2 = True
    else:
        isArray2 = False

#    print( "Telescopes %s,%s Max Dt = %.2f (s) " % \
#           (EventDir1['dir'], EventDir2['dir'], tOffset/OneMjdSec))
    # now match event times in directories 1 and 2
    for lll in range(n1):
        mjd1 = mjd1s[lll]
        # assume closest match is first file in 2nd list
        if isArray2: 
            t2, iii2 = findTime( mjd2s, mjd1)
        else:
            # else a single event, closest is that event
            t2 = float(mjd2s)
            iii2 = 0
        dt12 = abs( mjd1 - t2)
        if dt12 > tOffset:
            ii12s[lll] = NOMATCH
        else:
            ii12s[lll] = iii2
#            print(" %5d, match=%5d  dt=%8.2fs " % (lll, iii2, dt12*86400.))
        dt12s[lll] = dt12
        
    return ii12s, dt12s
# end of find pairs of matches

def showMatch ( matchIndex, aveMjd, nDir, matches, eventDirs):
    """
    showMatch() prints all matches for an event found
    where
    matchIndex == integer identifing this match
    nDir == number of telescopes
    matches == index to event matching in each telescope (length nDir)
    eventDirs == structure containing lists of events for each telescope
    """
    nMatch = 0
    for iDir in range( nDir):
        iMatch = matches[iDir]
        if iMatch != NOMATCH:
            nMatch = nMatch + 1
    if nMatch <= 0:
#        print("No Matches for Match Index: %5d" % (matchIndex))
        return
    print( "Match Index %5d (%d) Ave Mjd = %12.6f" % (matchIndex, nMatch, aveMjd))

    for iDir in range( nDir):
        iMatch = matches[iDir]
        if iMatch == NOMATCH:
            continue
        filei = eventDirs[iDir]['events'][iMatch]
        print( "file %2d: %5d: %s" % (iDir, iMatch, f"{filei:>45}"))
    return

def copyMatch ( matchIndex, outDir, aveMjd, nDir, matches, eventDirs):
    """
    copyMatch() copies flashes for a match to output directories
    where
    matchIndex == integer identifing this match
    nDir == number of telescopes
    matches == index to event matching in each telescope (length nDir)
    eventDirs == structure containing lists of events for each telescope
    """
    nMatch = 0
    for iDir in range( nDir):
        iMatch = matches[iDir]
        if iMatch != NOMATCH:
            nMatch = nMatch + 1
    if nMatch <= 0:
#        print("No Matches for Match Index: %5d" % (matchIndex))
        return
    print( "Match Index %5d (%d) Ave Mjd = %12.6f" % (matchIndex, nMatch, aveMjd))

    for iDir in range( nDir):
        iMatch = matches[iDir]
        if iMatch == NOMATCH:
            continue
        filei = eventDirs[iDir]['events'][iMatch]
        print( "copying %2d: %5d: %s" % (iDir, iMatch, f"{filei:>45}"))
        fileparts = filei.split('/')
        fileDir = outDir + "/" + fileparts[0]
        try:  # first try and create the directory
            os.makedirs(fileDir)
            print( "Creating directory: %s" % (fileDir))
        except:
            pass
        # now copy event flie
        shutil.copy(filei, fileDir)
    return

def showEvent ( nDir, event, eventDirs):

    """
    showEvent  prints the event info and match
    """
    
    matches = event['list']
    aveMjd = event['mjd']
    nEvent = event['n']
    showMatch( nEvent, aveMjd, nDir, matches, eventDirs)
    return

def showEvents ( nMatch, eventList, nDir, eventDirs):

    """
    showEvent  prints the event info and match
    """

    for iMatch in range( nMatch):
        event = eventList[iMatch]
        matches = event['list']
        aveMjd = event['mjd']
        showMatch( iMatch, aveMjd, nDir, matches, eventDirs)
    return

def listEvents ( iDir, nDir, eventDirs):

    """
    listEvents  prints a list of events
    """
    nEvent = eventDirs[iDir]['n']
    dirName = eventDirs[iDir]['dir']
    print("Dir: %2d nEvents: %5d in Directory: %s" % \
          (iDir, eventDirs[iDir]['n'], dirName))
    
    for iEvent in range (nEvent):
        mjd = eventDirs[iDir]['mjds'][iEvent]*86400.
        file = eventDirs[iDir]['events'][iEvent]
        print(" %5d %15.8f s: %s " % (iEvent, mjd, file ))
    return

def findMatches( nDir, eventDirs, telNames, tOffset = OneMjdSec, verbose = False):
    """
    findMatches() takes the list of events in telescopes=number of directories
    and creates a list of matches to events in other teleseops
    where
    nDir == number of telescopes active.
    eventDirs is a 3D structure of events found.  There is one element for each
    telescope the elements contain.  Each Element is two dimensional 
    1: number of events found
    2: name of event directory
    3: list of names of event files
    4: list of mjds of events found
    telNames == is a list of short names of telescopes
    The function returnes
    1: count of matches found
    2: list of event matches.
    """

    # Make a list of counts of events for each telescope
    nEves = np.zeros( nDir, dtype=int)

    print("      ", end="")
    # print labels and fill of tne nEvents list
    for iDir in range( nDir):
        nEves[iDir] = int(eventDirs[iDir]['n'])
        print(" %5s" % (telNames[iDir]), end = "")
    print("")
    # prepare to count matches
    listCount = 0
    grandTotal = 0
    nPair = 0
    ijmatrix=np.zeros((nDir,nDir),dtype=int)
    dtLists = []
    ijLists = []
    # now for all directories find the matches in other directories
    for iDir in range( nDir - 1):
        niEve = nEves[iDir]
        # for each of the pairs of directories, find closest matches
        mjds = eventDirs[iDir]['mjds']
        print("%5s " % (telNames[iDir]), end = "")

        countMatch = 0   # prepare to count the total number of matches to all other horns
        
        # now pad with blanks since a telescope does not match itself
        for iMatch in range( iDir+1):
            print( "      ", end = "")

        for jDir in range( iDir+1, nDir):
            njEve = nEves[jDir]
            # arrays iis and dts have length of the number of events in iDir
            ijs, dts = findpairs( eventDirs[iDir], eventDirs[jDir], tOffset=tOffset)
            # now create lists of pairs for eachs pair of telescopes.
            # This is a two dimentional list, with the iDir=0 to nDir-1
            # The second list is shorter iDir=1 to nDir,
            # if the first pass for telescope iDir
            ijmatrix[iDir][jDir] = listCount   # the index into the list is a complex function of i,j indices, just keep list
            ijLists.append( ijs)
            dtLists.append( dts)
            listCount = listCount + 1

            aCount = 0
            for iEve in range( niEve):
                if ijs[iEve] != NOMATCH:
                    aCount = aCount + 1
#            matchCounts[iDir] = aCount
            print( "%5d " % (aCount), end = "")
            countMatch = countMatch + aCount
        # updata total count of all matches for all horns
        grandTotal = grandTotal + countMatch
        print( "Total %6d" % (countMatch))
    ######################################################################

    # prepare to initial each event for each telescop to no match.
    # create an array with no matches in all directoreis
    noMatches = np.zeros(nDir, dtype=int) +  NOMATCH
    bigDts = np.zeros(nDir) + MAXDT
    zeroCounts = np.zeros(nDir, dtype=int)

    nMatch = 0
    eventList = [ ]
    
    # for all events in first directory, find acceptable matches in others
    for iDir in range(nDir - 1):
        niEve = nEves[iDir]                          # check each event found by this telescope
        mjdsi = eventDirs[iDir]['mjds']              # get list of all mjds for this telescope
        for iEve in range( niEve):
            matches = copy.deepcopy(noMatches)       # assume no matches to this event, except the telescope itself
            matchDts = copy.deepcopy(bigDts)        
            matches[iDir] = iEve                     # This is the event to be matched.
            matchDts[iDir] = 0.                      # no offset from the this event
        # matches is the length of the number of telescopes
            mjdi = mjdsi[iEve]
#            mjdi = eventDirs[iDir]['mjds'][iEve]      # get mjd for this event
#            matchCounts = copy.deepcopy( zeroCounts) # counting dublicate matches to event
#            matchCounts[iDir] = 1                    # this event only has one (exact) match
            for jDir in range( iDir+1, nDir):        # comparing with events in other telescopes
                ijIndex = ijmatrix[iDir][jDir]
                ijs = ijLists[ijIndex]               # ijs is a list of indexes of the 1st telescope, with matches to 2nd telescope
                ijdts = dtLists[ijIndex]             # time offsets for matches, often very large #s, never zero
                jMatch = int(ijs[iEve])              # j-event matchies this i-event for telescope iDir.
                if jMatch == NOMATCH:                # if no matches od ith event in the jDirectory,
                    continue                         # no need to test other events in j
                mjdj = eventDirs[jDir]['mjds'][jMatch]      # get mjd for this event
                dt = abs( mjdi - mjdj)               # time offset betwee iEve and jMatch
                if dt < tOffset:
                    matches[jDir] = jMatch
                    matchDts[jDir] = dt
                    nPair = nPair+1
                    if int(nPair/30000)*30000 == nPair and nPair > 0:
                        print( "Matching events: Telescopes %3d,%3d Events (%5d, %5d)" % \
                               (iDir, jDir, iEve, jMatch))
            # end for all other telescpes
            aveMjd = 0.
            nTel = 0
            for kDir in range( nDir):
                kMatch = matches[kDir]
                if kMatch != NOMATCH:
                    mjd = eventDirs[kDir]['mjds'][kMatch]
                    aveMjd = aveMjd + mjd
                    nTel = nTel + 1
            if nTel < 2:         # must match with at least one other telescope to be a match
                continue
            aveMjd = aveMjd/float(nTel)
            # event: 1: nmatch index number, 2: count of matching telescopes,
            #    3: average mjd and 4: list of matching event indicies
            notGroup = ' '
            event = { 'nmatch': nMatch, 'count': nTel, 'mjd': aveMjd, 'flag': notGroup, 'list': matches }
            if verbose:
                showMatch( nPair, aveMjd, nDir, matches, eventDirs)
            eventList.append( event)
            nMatch = nMatch + 1
        # end for all events
    #end for all telescopes
    if verbose:
        print( "Total number of event matches found %d" % (nMatch))
    return nMatch, eventList   # end of findMatches()

def sumMjds( matches, nDir, eventDirs, verbose = False):
    """
    sumMjds() recalculates the average mjd and nTel for a list of matches
    where
    matches is a list of nDir possible matches for an event
    nDir is the number of telescopes possible for this match
    eventDirs is the main data structure for all events in the group
    """
    nTel = 0
    aveMjd = 0.
    for iDir in range( nDir):
        iMatch = matches[iDir]
        if iMatch != NOMATCH:
            mjd = eventDirs[iDir]['mjds'][iMatch]
            aveMjd = aveMjd + mjd
            nTel = nTel + 1
    if nTel > 0:
        aveMjd = aveMjd / float( nTel)
    return nTel, aveMjd

def trimEvents( nEvents, events, nDir, eventDirs, verbose = False):
    """
    trimEvents: read through all matched events and removes duplicate matches of a telescope event in
    in multiple matched events.  This puts an event in a telescope in only one match.
    This clears up the match of a single event in one telescope to many events in another telescope
    """
    # Make a list of counts of events for each telescope
    nEves = np.zeros( nDir, dtype=int)

    outEvents = copy.deepcopy( events)
    for iDir in range( nDir):
        # now get the mjd of this event for this telescope
        nEves[iDir] = int(eventDirs[iDir]['n'])

    for iEvent in range( nEvents-1):
        eventi = events[iEvent]
        matchesi = eventi[ 'list' ]
        neventi  = eventi['nmatch']
        aveMjdi = eventi[ 'mjd' ]
        for iDir in range( nDir):
            iMatch = matchesi[iDir]
            if iMatch == NOMATCH:
                continue
            # now get the mjd of this event for this telescope
            eventMjd = eventDirs[iDir]['mjds'][iMatch]
            dti = abs( eventMjd - aveMjdi)
            ni = eventi['count']
            for jEvent in range (iEvent+1, nEvents):
                eventj = events[jEvent]
                matchesj = eventj['list']
                jMatch = matchesj[iDir]
                if iMatch != jMatch:   # not the same matching event in telescope i
                    continue
                nj = eventj['count']
                aveMjdj = eventj['mjd']
                dtj = abs( eventMjd - aveMjdj)
                if verbose:
                    print("Checking Directory %d matches %d and %d (%d, %d)" % \
                          (iDir, iEvent, jEvent, ni, nj))
                # only keep one of the two matches, flag the other as NOMATCH
                if ni == nj:         # if the same number of matches
                    if dtj < dti:    # keep the closest event to the average
                        keepj = True
                    else:
                        keepj = False
                elif nj > ni:        # else if more matches assume that is the best match
                    keepj = True
                else:
                    keepj = False
                # if keeping the jth event
                if keepj:
                    matchesi[iDir] = NOMATCH
                    nTel, aveMjd = sumMjds( matchesi, nDir, eventDirs)
                    if verbose:
                        print( "Ni matches goes from %d to %d (%d,%d)" % \
                               (ni, nTel,neventi, neventj))
                    eventi['list'] = matchesi
                    eventi['mjd']  = aveMjd
                    eventi['count'] = nTel
                    outEvents[ iEvent] = eventi
                else:  # flagging j instead of keeping
                    matchesj[iDir] = NOMATCH
                    nTel, aveMjd = sumMjds( matchesj, nDir, eventDirs)
                    if verbose: 
                        print( "Nj matches goes from %d to %d (%d,%d)" % \
                               (nj, nTel, neventi, neventj))
                    eventj['list'] = matchesj
                    eventj['mjd']  = aveMjd
                    eventj['count'] = nTel
                    outEvents[ jEvent] = eventj
                    # if replacing the jEvent, then closer matches are still possible
            # end for all other events
        # end for all events
    # end for all directories
    # return the newly trimmed events
    return outEvents
    
def compressEvents( nEvents, events, nMin = 3, verbose = False):
    """
    compressEvents() takes a list of events and removes the matches with too few telscopes
    """
    nOut = 0
    eventsOut = []
    for iEvent in range( nEvents):
        event = copy.deepcopy(events[iEvent])
        nTel = event['count']
        if nTel >= nMin:
            eventsOut.append( event)
            nOut = nOut + 1

    if verbose:
        print( "%5d input events trimed to %5d output events (delta %d)" % \
               (nEvents, nOut, nEvents-nOut))
    return nOut, eventsOut
    # end of compress Events
    
if __name__ == "__main__":
    """
    Create a set of events and find matches
    """
    import random
    tOffset = OneMjdSec
    mjd0 = 0.    # default mjd

    telNames = ["pi1", "pi2", "pi3", "pi4", "pi5"]
    nDir = len (telNames)
    print( "Creating simulated list of Events for %d telescopes" % (nDir))
    print( "Each Telescope finds a different number of events")
    print( "With different event times.  Match offset: %.2f (s)" % (tOffset*86400.))
    

    for iDir in range( nDir):
        n = iDir + 5
        mjds = np.zeros(n) + mjd0                                     
        for iName in range( n):
            mjds[iName] = (float(iName+1.)*(OneMjdSec*3.)) + random.gauss(0., (iName+1.)*.2)*OneMjdSec
        eventNames = []
        for iName in range( n):
            eventNames.append(  str(telNames[iDir]+("_%d_%.3f" % (iName, mjds[iName]*86400.))))
            eventDir =  { 'dir': telNames[iDir], # name of directory
                      'events' : eventNames,     # names of events in directories
                      'n': n,                    # count of events in directory
                      'mjds': mjds,              # array of MJDs of these events
                     }
        if iDir == 0:
            eventDirs = { iDir: eventDir }
        else:
            eventDirs.update( { iDir: eventDir })

    # confirm events are initialized correctly
    listEvents( 1, nDir, eventDirs)

    nMatch, eventList = findMatches( nDir, eventDirs, telNames, tOffset = tOffset, verbose = True)
    # End of test module

    print("##############################################################################")
    print("Starting Event triming to remove duplicate matches to a single telescope-event")
    
    eventTrim = trimEvents( nMatch, eventList, nDir, eventDirs, verbose = False )

    nMin = 3
    nRemain, eventFinal = compressEvents( nMatch, eventTrim, nMin=nMin, verbose = True)
    print("##############################################################################")
    print("Compressed %d events to %d events with minimum number of telescopes=%d" % \
          (nMatch, nRemain, nMin))
    
    # now show the remaining events
    showEvents( nRemain, eventFinal, nDir, eventDirs)
