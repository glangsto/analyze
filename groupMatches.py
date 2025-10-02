#Python find group matches into categories 
#HISTORY
#25Oct02 GIL add gaussian fitting
#25Sep29 GIL initial version from code in matchevents.py

import sys
import os
import copy
import numpy as np
from findMatches import *

def groupMatches( mjdRef, nMatch, nDir, eventLists, nDay, verbose = False):
    """
    groupMatches takes a set of input events and categorizes them based on counts of events in an interval
    where
    mjdRef - Reference Date for matching.
    nMatch - Total number of events in the input list
    nDir   - Number of telescopes with matches
    eventLists - list of events each with multiple matches to up to nDir telescopes
    nDay   - Number of parts of a day for grouping matches.
    verbose - optionally flag printing diagnostic info.
    """

# event: 1: nmatch index number, 2: count of matching telescopes,                                                       
#    3: average mjd and 4: list of matching event indicies                                                              
    notGroup = " "
# event = { 'nmatch': nMatch, 'count': nTel, 'mjd': aveMjd, 'flag': notGroup, list': matches }

    # show first and last matches 
    if verbose:
        matches = eventLists[0]['list']
        aveMjd = eventLists[0]['mjd']
        dMjd = aveMjd - mjdRef
        iMatch = eventLists[0]['nmatch']
        print( "%5d %6.2f" % (iMatch, dMjd))
        print( matches)

    # prepare to count events in intervals
    countInDt = np.zeros(nDay, dtype=int)
    # create an empty group for each fraction of a day
    groupOfGroups = [[] for _ in range(nDay)]
    
    # for all matches, is this event in this time range
    for iMatch in range(nMatch):
        mjd = eventLists[iMatch]['mjd']
        dMjd = mjd - mjdRef
        # compute index into counts as a fraction of a day.
        iDay = int((dMjd*nDay))
        # wrap around day fraction
        iDay = iDay % nDay
        # add this index to the lost of groups for this time range
        groupOfGroups[iDay].append( iMatch)
        countInDt[iDay] = countInDt[iDay] + 1
        
    iGroup = 0             # initialize count of groups
    groupLists = [ ]       # initialize output lists of event groups
    nTen = 0
    nHundred = 0
    nThousand = 0
    dMjd = 1./nDay
    # for each fraction of a day
    for iDay in range( nDay):
        # if any events in this fraction
        if countInDt[iDay] > 0:
            iFirstEvent = groupOfGroups[iDay][0]
            if verbose:
                print("%7.2f (%3d) %5d first Event=%5d" % (iDay*dMjd, iDay, countInDt[iDay], iFirstEvent))
            # index to first match in list for this group
            # initialize output group with 1st event in group
            aveMjd = 0.
            maxEvent = iFirstEvent
            for iEvent in range (countInDt[iDay]):
                inMatch = groupOfGroups[iDay][iEvent]
                # keep the event with most matches to represent the group
                if eventLists[inMatch]['count'] > eventLists[maxEvent]['count']:
                    if verbose:
                        print( "Event %5d(%5d) has %3d Telescopes > Event %5d with %3d Telescopes" % \
                               (inMatch, iEvent, eventLists[inMatch]['count'], \
                               maxEvent, eventLists[maxEvent]['count']))
                    maxEvent = inMatch
                # still averaging all times for all events in range
                aveMjd = aveMjd + eventLists[inMatch]['mjd']
            # use the maximum telescoope match event for comparision of types of events in group
            nTel = eventLists[maxEvent]['count']
            groupEvent = copy.deepcopy(eventLists[maxEvent]) 
            matches = groupEvent['list']
            count = 0
            # make sure this representative event has some matches
            for iDir in range( nDir):
                if matches[iDir] != NOMATCH:
                    count = count + 1
            if count < 2:
                print("Few matches for this representative Event: %d %d %d" % \
                      (iEvent, maxEvent, count))
                continue
                      
            if verbose and nTel > 0:
                print("Max telescopes %d in time %3d, event %d" %
                      (nTel, iDay, maxEvent))
            groupEvent['count'] = countInDt[iDay]
            groupEvent['nmatch'] = iGroup
            groupEvent['mjd'] = aveMjd/countInDt[iDay]
            if countInDt[iDay] == 1:
                groupEvent['flag'] = " "
            elif countInDt[iDay] < 10:
                groupEvent['flag'] = ("%d" % countInDt[iDay])
            elif countInDt[iDay] < 100:    # else use Roman numerals for catagories
                groupEvent['flag'] = "X"
                nTen = nTen + 1
            elif countInDt[iDay] < 1000:
                groupEvent['flag'] = "C"
                nHundred = nHundred + 1
            else:
                groupEvent['flag'] = "M"
                nThousand = nThousand + 1
            iGroup = iGroup + 1
            
            groupLists.append( groupEvent)
        # end if any events in this fraction of a day
    # end for all fractions of a day

    return iGroup, nTen, nHundred, nThousand, groupLists
# end of groupMatches

def main( ):
    """
    Test function
    """

    nDir = 20
    nMatch = 1000
    eventLists = []
    iMatch = 0
    notGroup = ' '
    matches = np.zeros( nDir, dtype=int)
    aveMjd = 0.
    nDay = 144
    dT = 1./nDay
    for iDir in range(nDir):
        t = 0.
        matches = np.zeros( nDir, dtype=int)
        for iMatch in range (nMatch):
            matches[iDir] = 1
            t = iMatch*dT/3. + np.random.normal(0.,.05)
            event = { 'nmatch': iMatch, 'count': iDir, 'mjd': t, 'flag': notGroup, 'list': matches }
            eventLists.append( event)

    nGroups, nTen, nHundred, nThousand,groupLists = groupMatches( aveMjd, nMatch, nDir, eventLists, nDay, verbose = True)

    print( "Matched events into %d groups" % (nGroups))
    print( "%d groups with > 10 matches, %d groups with > 100 matches and %d with > 1000" % \
           (nTen, nHundred, nThousand))
    return
                  
if __name__ == "__main__":
    main()
    
