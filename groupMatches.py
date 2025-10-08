#Python find group matches into categories 
#HISTORY
#25Oct05 GIL clean up fit example
#25Oct03 GIL create new new list of fits to groups
#25Oct02 GIL add gaussian fitting
#25Sep29 GIL initial version from code in matchevents.py

import sys
import os
import copy
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
from findMatches import *

def gaussian(x, amplitude, mean, std_dev):
    return amplitude * norm.pdf(x, mean, std_dev)
   
def groupFit( nDay, inData, nGroup, groupList, nFitMax = 3, verbose = False):
    """
    groupFit takes a list of groups and fits a gaussian to the time verses number of events.
    where
    nDay      = number of fractions of a day for counts
    inData    = input array of events as a function of time
    nGroup    = number of groups in list, will augment coordiate times
    groupList = list of groups to be fit
    nFitMax   = maximum number of fits to attempt
    returns
    nFit - number
    fits - a three component set of fit values 
        0: time of center of fit,
        1: peak of group count
        2: FWHM of the fit time range for fit.  0 indicates no Fit
    """

    aFit = { 'time': 0.0, 'peak': 0.0, 'stdDev': 0.0, 'fwhm': 0.0, 'rms': 0.0, 'goodFit': False }
    noFit = copy.deepcopy( aFit)

    fits = [ ]
    if nDay < 24:
        print( "Number of intervals (%d) too small for a fit to event counts" % (nDay))
        nFit = 0
        fits.append( noFit)
        return nFit, fits
    
    # prepare to set the time axis
    dt = 1./nDay   # set the time interval for each box
    times = np.zeros(nDay)
    
    iMax = 0
    time = dt/2.   # set middle of time range
    # initialize the mjd time array
    for iDay in range (nDay):
        times[iDay] = time
        time = time + dt

    # fill the scratch dataset
    data = copy.deepcopy( inData)

    outFit = 0
    # now for up to the maximun number of fits, t
    for iFit in range(nFitMax):
        iMax = 0
        for iDay in range( nDay):
            if data[iDay] > data[iMax]:
                iMax = iDay
            
        # Initial guesses for parameters (optional but recommended for better fits)
        initial_guesses = [data[iMax], times[iMax], np.std( times)]
        if verbose:
            print( "Fitting Gaussian to %d groups, max %d at time %f h (%d)" % \
                   (nGroup, data[iMax], times[iMax]*24., iMax))

        # Fit the Gaussian
        try:
            #    params, covariance = curve_fit(gaussian, bin_centers, hist_counts, p0=initial_guesses)
            # initial_guesses = [max(bin_heights), np.mean(data), np.std(data)]
            # popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, p0=initial_guesses)
            # Extract fitted parameters
            # namplitude_fit, mean_fit, std_dev_fit = popt

            params, covariance = curve_fit(gaussian, times, data, p0=initial_guesses)
        except:
            break
        amplitude, meanTime, stdDev = params
        aFit = { 'time': 0.0, 'peak': 0.0, 'stdDev': 0.0, 'fwhm': 0.0, 'rms': 0.0, 'goodFit': False }
        # 5. Calculate FWHM
        fwhm = 2 * np.sqrt(2 * np.log(2)) * stdDev
        aFit['time'] = meanTime
        aFit['fwhm'] = fwhm
        aFit['stdDev'] = stdDev
        aFit['goodFit'] = True
        # the fit returns parameters not directly related to peak amplitude
        y_fit = gaussian(times, *params)
        # use the calculated values to find peak amplitude
        amplitude = np.max(y_fit)
        aFit['peak'] = amplitude
        # now remove the fit just found, so the next peak can be found
        for iDay in range( nDay):
            if verbose and data[iDay] > 1.:
                print("%3d %12.2f %12.2f %12.2f" % (iDay, data[iDay], y_fit[iDay],  data[iDay] - y_fit[iDay]))
            data[iDay] = data[iDay] - y_fit[iDay]
        # use the residuals to compute estimated RMS of a fit.  Includes noise in remainging events.
        rms = np.sqrt(np.mean(np.square(data)))
        aFit['rms'] = rms
        if verbose:
            print( "Fit: peak: %.2f+/-%.2f, at Mjd: %f +/- %f (h)" % \
                   (amplitude, rms, meanTime*24, stdDev*24))
        # if peak is significant compared to RMS
        if amplitude > 5. * rms:
            # here a good fit
            fits.append( aFit)
            outFit = outFit + 1
        else:
            # else not a significant fit
            break
    # end for all fit attempts

    if outFit < 1:
        fits.append( noFit)
        print( "Gausian fitting Failed!")
            
    return outFit, fits
   # end of groupFit
            
def groupMatches( mjdRef, nMatch, nDir, eventLists, nDay, verbose = False):
    """
    groupMatches takes a set of input events and categorizes them based on
    counts of events in an interval
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

    if nMatch < 1:
        print("No Mathes to group")
        groups = []
        return 0, 0, 0, 0, groups

    # show first and last matches 
    if verbose:
        matches = eventLists[0]['list']
        aveMjd = eventLists[0]['mjd']
        dMjd = aveMjd - mjdRef
        iMatch = eventLists[0]['nmatch']
        print( "%5d %6.2f" % (iMatch, dMjd))
        print( matches)

    # prepare to count events in intervals
    countsInDt = np.zeros(nDay, dtype=int)
    times = np.zeros(nDay)      # time of each event interval
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
        countsInDt[iDay] = countsInDt[iDay] + 1
        
    iGroup = 0             # initialize count of groups
    groupLists = [ ]       # initialize output lists of event groups
    nTen = 0
    nHundred = 0
    nThousand = 0
    dMjd = 1./nDay
    # for each fraction of a day
    for iDay in range( nDay):
        # if any events in this fraction
        if countsInDt[iDay] > 0:
            iFirstEvent = groupOfGroups[iDay][0]
            if verbose:
                print("%7.3f (%3d) %5d first Event=%5d" % (iDay*dMjd, iDay, countsInDt[iDay], iFirstEvent))
            # index to first match in list for this group
            # initialize output group with 1st event in group
            aveMjd = 0.
            maxEvent = iFirstEvent
            for iEvent in range (countsInDt[iDay]):
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
            groupEvent['count'] = countsInDt[iDay]
            groupEvent['nmatch'] = iGroup
            groupEvent['mjd'] = aveMjd/countsInDt[iDay]
            if countsInDt[iDay] == 1:
                groupEvent['flag'] = " "
            elif countsInDt[iDay] < 10:
                groupEvent['flag'] = ("%d" % countsInDt[iDay])
            elif countsInDt[iDay] < 100:    # else use Roman numerals for catagories
                groupEvent['flag'] = "X"
                nTen = nTen + 1
            elif countsInDt[iDay] < 1000:
                groupEvent['flag'] = "C"
                nHundred = nHundred + 1
            else:
                groupEvent['flag'] = "M"
                nThousand = nThousand + 1
            iGroup = iGroup + 1
            
            groupLists.append( groupEvent)
        # end if any events in this fraction of a day
    # end for all fractions of a day
    nGroup = iGroup
    

    return nGroup, nTen, nHundred, nThousand, groupLists
# end of groupMatches

def main( ):
    """
    Test function
    """

    import random
    import matplotlib.pyplot as plt

    # assume 30 telescopes
    nDir = 30
    nMatch = 10000
    eventLists = []
    iMatch = 0
    notGroup = ' '
    matches = np.zeros( nDir, dtype=int)
    aveMjd = 0.
    nDay = 144
    dT = 1./nDay
    totalCountsInDt = np.zeros(nDay, dtype=int)
    times = np.zeros(nDay)

    dt = 1./nDay
    t = dt/2.
    for iDay in range(nDay):
        times[iDay] = t
        t = t + dt
        
    tEvent = .5   # noon mjd
    fwhm  = 20./1440.   # +/- 20 minute sof a day
    # for all simulated telescopes
    random.seed(42)  # for reproducibility
    for iDir in range(nDir):
        t = 0.
        matches = np.zeros( nDir, dtype=int)
        # for all simulated matches
        for iMatch in range (nMatch):
            matches[iDir] = 1
            t = random.gauss( tEvent, fwhm)
            event = { 'nmatch': iMatch, 'count': iDir, 'mjd': t, 'flag': notGroup, 'list': matches }
            eventLists.append( event)
            # compute index into counts as a fraction of a day.
            iDay = int(t*nDay)
            # wrap around day fraction
            iDay = iDay % nDay
            # now count events in this time interval
            totalCountsInDt[iDay] = totalCountsInDt[iDay] + 1
                            
    nGroups, nTen, nHundred, nThousand,groupLists = groupMatches( aveMjd, nMatch, nDir, eventLists, \
                                                                  nDay, verbose = True)

    print( "Matched events into %d groups" % (nGroups))
    print( "%d groups with > 10 matches, %d groups with > 100 matches and %d with > 1000" % \
           (nTen, nHundred, nThousand))

    # now check if any groups are suitable for fitting
    if nHundred > 0 or nThousand > 0:
        nFit = 2
        nFit, fits = groupFit( nDay, totalCountsInDt, nGroups, groupLists, nFitMax=nFit, verbose = True)
        if nFit < 1:
            print("No Fits found to event counts")
            return
        else:
            print("%d Fits found to event counts" % (nFit))

    popt = np.zeros(3)
    popt[0] = fits[0]['peak']
    popt[1] = fits[0]['time']
    popt[2] = fits[0]['fwhm']
    stdDev  = fits[0]['fwhm']
    mean = fits[0]['time']
    fwhm = fwhm = 2 * np.sqrt(2 * np.log(2)) * stdDev
             
    y_fit = gaussian(times, *popt)
    # for plotting convert times to hours
    times = times*24.
    mean = mean * 24.
    stdDev = stdDev * 24.
    fwhm = fwhm * 24.
    plt.plot(times, totalCountsInDt, 'g.', linewidth=9, markersize=12, label=f'Input Counts')
    plt.plot(times, y_fit, 'r-', label=f'Gaussian Fit:\nMean={mean:.2f}, Std Dev={stdDev:.2f}\nFWHM={fwhm:.2f}')
    plt.xlim( 10., 14.)
    
    # 7. Add plot details
    plt.title('Gaussian Fit to Histogram with FWHM')
    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"Fitted Mean: {mean:.2f} hours")
    print(f"Fitted Standard Deviation: {stdDev:.2f}")
    print(f"Calculated FWHM: {fwhm:.2f} hours")
    return
                  
if __name__ == "__main__":
    main()
    
