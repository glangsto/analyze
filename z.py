##Python Script to find az and el offsets based on galactic plan zero latitude crossings
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#16Oct10 GIL simplify
#16Aug03 GIL improve iterated search
#16Aug02 GIL test for finding az,el offsets
#16Oct07 GIL check for changes in frequency and/or gain
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import radioastronomy
import copy
import interpolate
import Moment
import os.path

avetimesec = 3600.
dy = -1.
linelist = [1420.0, 1419.0, 1418.0]  # RFI lines in MHz
linewidth = [11, 5, 9]
firstrun = True

nargs = len(sys.argv)

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1e8

c = 299792.  # (v km/sec)
nuh1 = 1420.4056 # neutral hydrogen frequency (MHz)
thot = 285.0  # define hot and cold
#thot = 272.0  # 30 Farenheit = 272 K
tcold = 10.0
tmin = 20.0 
tmax = 999.0 # define reasoanable value limits

# prepare to remove a linear baseline
xa = 550
#xb = 1024-xa
xb = 1024-200

nplot = 0

#first argument is the averaging time in seconds

avetimesec = float(sys.argv[1])
print "Average time: ", avetimesec, " (seconds)"

# first read through all data and find hot load
names = sys.argv[2:]
names = sorted(names)
minvel = -200.
maxvel = 200.
minvel = -130.
maxvel = 130.

minvel = -150.
maxvel = 150.

eloffset = -.77
azoffset = -34.04

nprint = 0

def interpolate_range( minvel, maxvel, vel, ys):
    """
    Interpolate arrays over a range of velocities
    minvel and maxvel are the ranges to interpolate
    input vel array of velocities 
    input ys values to interpolate
    """

    imin, imax = velocity_to_indicies( vel, minvel, maxvel)
        
    ya = np.median(ys[(imin-10):(imin+10)])
    yb = np.median(ys[(imax-10):(imax+10)])
    slope = (yb-ya)/(imax-imin)
    # finally interpolate over values
    for iii in range( imin, imax):
        ys[iii] = (ya + (slope*(iii-imin)))
    
    return ys

def compute_sum( minvel, maxvel, ra, dec, gallon, gallat, vel, tsky, n):
    """ 
    compute the integral of an array over a relevant velocity range
    """
    global firstrun
    nvel = len(vel)
    nsky = len(tsky)
    if (nsky > nvel):
        print "Error writing sum: Velocites > Intensity array sizes: ", nvel, nsky

    if nvel < 5:
        print 'Not enough values in arrays to compute reliable sums'

    imin, imax = velocity_to_indicies( vel, minvel, maxvel)  

    if firstrun:
        print 'Velocity range: %9.2f to %9.2f km/sec (%d,%d)' % (minvel, maxvel,imin,imax)
        firstrun = False

# first get the RMS in a region not expected to have signal
    vs = vel[imin-20:imin]
    ts = tsky[imin-20:imin]
    (ave,adev,rms,var,skew,curt) = Moment.moment(vs, ts)    

# now get the sum over the range of interest
    vs = vel[imin:imax]
    ts = tsky[imin:imax]
    (ave,adev,sdev,var,skew,curt) = Moment.moment(vs, ts)    

# now compute the intensity weighted velocity sum.
    vts = vel[imin:imax]*tsky[imin:imax]

    # compute the weighted average velocity
    (vave,vadev,vsdev,vvar,vskew,vcurt) = Moment.moment(vs, vts)    
    # normalize for average temperature
    if ave == 0: 
        print 'Average temperature is zero, can not normalize average velocity'
    else:        
        vave = vave/ave
        vsdev = vave/ave

# File = 2017-01-03.sum
# RA  DEC   T_sum  T_std   V_ave  V_std   Time  
# 54.56  -35.32    -5.587     3.039  -624.929   111.851 20:02:38
    tsum = 0
    for iii in range(imax-imin):
        tsum = tsum+ts[iii]
# take absolute value of velocity
    iref = int(nvel/2)
    vref = vel[iref]
    dv   = (vel[iref+2]-vel[iref-2])/4.
    if dv < 0.:
        dv = -dv
    tsum = tsum * dv  # scale integral by velocity resolution
    return tsum

def findhotgain( names, avetimesec):
    lastgain = 0.
    # flag start of search for a latitude zero crossing
    ncrossings = -1
    crossingsum = 0.
    lastgallat = -100.
    lastsum = 0.

    minGlat = +90.
    maxGlat = -90.
    nhot = 0         # number of obs with el < 0
    ncold = 0
    nhigh = 0        # highest galactic latitude
    lastfreq = 0.
    lastbw = 0.
    lastgain = 0.
    lastel = 0.
    lastaz = 0.
    firstrun = True

    for filename in names:

        rs = radioastronomy.Spectrum()
        rs.read_spec_ast(filename)
        rs.telel = rs.telel + eloffset
        rs.telaz = rs.telaz + azoffset
        rs.azel2radec()    # compute ra,dec from az,el

        if rs.gains[0] != lastgain:
            if lastgain != 0:
                print 'Gain Change: ', lastgain, ' -> ', rs.gains[0]
                lastgain = rs.gains[0]

        if rs.telel < 0:
            if nhot == 0:
                hot = copy.deepcopy( rs)
                nhot = 1
            else:
                hot.ydataA = hot.ydataA + rs.ydataA
                hot.count = hot.count + rs.count
                nhot = nhot + 1
        else: # else above horizon, find min, max galactic latitudes
            if rs.gallat > maxGlat:
                maxGlat = rs.gallat
            if rs.gallat < minGlat:
                minGlat = rs.gallat

    if nhot > 0:
        hot.ydataA = scalefactor * hot.ydataA / float(nhot)
        print "Found %3d hot load obs" % nhot
    else:
        print "No hot load data, can not calibrate"
        exit()

    xv = hot.xdata * 1.E-6
    yv = hot.ydataA
    nData = len(xv)
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

    print 'Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % (minGlat, maxGlat)

# assume only a limited range of galactic latitudes are available
# not range about +/-60.
    use60Range = False

# all galactic latitudes above +/-60d can be used
    if minGlat < -60. or maxGlat > 60.:
        minGlat = -60.
        maxGlat = 60.
    else: # else no high galactic latitude data
    # use highest galactic latitudes - +/-5.degrees
        if -minGlat > maxGlat:  # if negative latitudes higher
            minGlat = minGlat + 5.
            maxGlat = 90.
        else: # else positive latitudes higher
            maxGlat = maxGlat - 5.
            minGlat = -90.

    outnames = names
    nout = 0
# now average coldest data for calibration
    for filename in names:

        rs = radioastronomy.Spectrum()
        rs.read_spec_ast(filename)
        rs.telaz = rs.telaz + azoffset
        rs.telel = rs.telel + eloffset
        rs.azel2radec()    # compute ra,dec from az,el

        if rs.telel < 0:
            continue
        if rs.gallat > maxGlat or rs.gallat < minGlat:
            if nhigh == 0:
                high = copy.deepcopy( rs)
                high.ydataA = rs.ydataA
                nhigh = 1
            else:
                high.ydataA = high.ydataA + rs.ydataA
                high.count = high.count + rs.count
                nhigh = nhigh + 1
        else:
            # if here, then only these names will be used to find zero crossings
            outnames[nout] = filename
            nout = nout + 1

    if nhigh < 1.:
        print "No high galactic latitude data: can not calibrate"
        exit()
    else:
        high.ydataA = scalefactor * high.ydataA/nhigh
        print "Found %d High Galactic Latidue spectra" % (nhigh)
        yv = high.ydataA
        cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

# finally compute gain on a channel by chanel basis
    gain = np.zeros(nData)
    vel = np.zeros(nData)
    for iii in range(nData):
        gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)
        vel[iii] = c * (nuh1 - xv[iii])/nuh1

# now interpolate over galactic velocities
    gain = interpolate_range( minvel, maxvel, vel, gain)

    if nout > 0:
        outnames = outnames[0:(nout-1)]
    else:
        print 'No low latitude file names left!'
        exit()
    return gain, vel, outnames

def velocity_to_indicies( vel, minvel, maxvel):
    """
    Function to compute indices from velocity array and target velocities
    """

    nData = len(vel)
    iref = int(nData/2)
    
    vref = vel[iref]
    dv   = (vel[iref+2]-vel[iref-2])/4.
    imin = int(((minvel - vref)/dv) + iref)
    imax = int(((maxvel - vref)/dv) + iref) + 1

    if imax < imin: # swap indices if frequency opposit velocity
        temp = imin
        imin = imax
        imax = temp

    if imin < 0:
        print 'Imin Error computing baseline: ', imin
        imin = 0
    if imin >= nData:
        print 'Imin Error computing baseline: ', imin
        imin = nData-1

    if imax < 0:
        print 'Imax Error computing baseline: ', imax
        imax = 0
    if imax >= nData:
        print 'Imax Error computing baseline: ', imax
        imax = nData-1

# some channels will be selected using x[imin:imax] so,
# now make sure range increases
    if imin > imax:
        temp = imin
        imin = imax
        imax = temp
    return imin, imax

def fit_baseline( xs, ys, imin, imax, nchan):
    """
    fit baseline does a polynomical fit over channels in a select range
    The baseline is returned.
    """
    
    xfit = np.concatenate( (xs[imin-nchan:imin+nchan],xs[imax-nchan:imax+nchan]))
    yfit = np.concatenate( (ys[imin-nchan:imin+nchan],ys[imax-nchan:imax+nchan]))
# calculate polynomial
    z = np.polyfit(xfit, yfit, 3)
    f = np.poly1d(z)

# calculate new x's and y's
    yout = f(xs)

    return yout

def findzero( names, gain, vel, azoffset, eloffset, avetimesec):
    lastgain = 0.
    global nprint
    # flag start of search for a latitude zero crossing
    ncrossings = -1
    crossingsum = 0.
    lastgallat = -100.
    lastsum = 0.

    nData = len(vel)
    minGlat = +90.
    maxGlat = -90.
    nhot = 0         # number of obs with el < 0
    ncold = 0
    nhigh = 0        # highest galactic latitude
    lastfreq = 0.
    lastbw = 0.
    lastgain = 0.
    lastel = 0.
    lastaz = 0.
    firstrun = True

    avetime = datetime.timedelta(seconds=avetimesec)
    firstgallon = -1.
    ncrossings = 0

    nread = 0        
# now read through all data and average cold sky obs
    for filename in names:

        parts = filename.split('/')
        nparts = len(parts)
        aname = parts[nparts-1]
        parts = aname.split('.')
        aname = parts[0]
        extension = parts[1]
        parts = aname.split('T')
        date = parts[0]
        time = parts[1]
        time = time.replace('_', ':')  # put time back in normal hh:mm:ss format

# exclude hot load data for averaging
        if extension == 'hot':
            continue

        rs = radioastronomy.Spectrum()
#  print filename
        rs.read_spec_ast(filename)
        rs.telaz = rs.telaz + azoffset
        rs.telel = rs.telel + eloffset
        rs.azel2radec()    # compute ra,dec from az,el

# if a sky observation and close to galactic plane
        if rs.telel > 0. and (rs.gallat > -5.0 and rs.gallat < 5.0) :

# if first time reading data, set obs parameters
            if lastfreq == 0.:
                lastfreq = rs.centerFreqHz 
                lastbw = rs.bandwidthHz
                lastgain = rs.gains[0]
                lastaz = rs.telaz
                lastel = rs.telel
                cold = copy.deepcopy( rs)
                ncold = 0
                timesum = 0.
                firstutc = cold.utc
                lastutc = cold.utc

            newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
            if newAzEl and nprint < 5:
                print 'Must keep az,el constant to find pointing offsets!'
                print 'Last %7.2f,%7.2f; New: %7.2f,%7.2f' % (lastaz, lastel, rs.telaz, rs.telel)
                lastaz = rs.telaz
                lastel = rs.telel
                nprint = nprint + 1
                break

            if ncold > 1:
            # time difference is between mid-points of integrations
                dt = rs.utc - cold.utc 
            # add the time since midpoint of latests
                dt = dt + datetime.timedelta(seconds=rs.durationSec)
                lastutc = rs.utc

           # if time to average (or end of all files)
                if (dt > avetime) or (filename == sys.argv[nargs-1]):
                    cold.ydataA = cold.ydataA/float(timesum)
                # have complicated steps to simple get the average time.
                    deltatime = endtime - starttime
                    aveutc,duration = radioastronomy.aveutcs( firstutc, lastutc)
                    cold.utc = aveutc
                    cold.azel2radec() # recompute coordinates for average time
                    ra = cold.ra
                    dec = cold.dec
                    gallat = cold.gallat
                    gallon = cold.gallon

                    if ncrossings < 0:
                        lastgallat = gallat
                        ncrossings = 0
                        crossingsum = 0.
                        
                    xv = cold.xdata * 1.E-6
                    yv = cold.ydataA * scalefactor
                    yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
                    xmin = min(xv)
                    xmax = max(xv)
                    count = cold.count
                    note = cold.noteA
                    #print('%s' % note)
                    ncolor = min(nmax-1, nplot) 
                
                    tsky  = np.zeros(nData)    # initialize arrays
                    for jjj in range (0, nData):
                        tsky[jjj] = yv[jjj]/gain[jjj]
                        
                    imin, imax = velocity_to_indicies( vel, minvel, maxvel)
                    
                    ymed = np.median(tsky[imin:imax])
                    ya = np.median(tsky[(imin-10):(imin+10)])
                    yb = np.median(tsky[(imax-10):(imax+10)])
                    slope = (yb-ya)/(imax-imin)
                    nFit = 20
                    baseline = fit_baseline( vel, tsky, imin, imax, nFit)
                    tsky = tsky - baseline

                    ymin = min(tsky[imin:imax])
                    ymax = max(tsky[imin:imax])
                    thesum = compute_sum( minvel, maxvel, ra, dec, gallon, gallat, vel, tsky, ncold)
                    
                # finallly if this is a latitude crossing sum.
                    if ((lastgallat < 0. and gallat > 0.) or (lastgallat > 0. and gallat < 0.)) and lastgallat > -100.:
                        if abs(gallat) < 5.:
                    # if first crossing, init sum
                            if ncrossings == 0:
                                ncrossings = 1
                                crossingsum = thesum 
#                                print 'Zero Latitude: %7.1f %7.1f: %10.0f, %10.0f' % (lastgallat, gallat, lastsum, thesum)
                            else:
                                # else average crossings
                                ncrossings = ncrossings + 1
                                crossingsum = crossingsum + thesum 
                        if abs(lastgallat) < 5.:
                            if ncrossings == 0:
                                ncrossings = 1
                                crossingsum = lastsum
#                                print 'Zero Latitude: %7.1f %7.1f: %10.0f, %10.0f' % (lastgallat, gallat, lastsum, thesum)
                            else:
                                # else average crossings
                                ncrossings = ncrossings + 1
                                crossingsum = crossingsum + lastsum 

                    lastgallat = gallat
                    lastsum = thesum
                    # computation done, start sum over gain
                    ncold = 0

# if this was a new obs; restart the sums
            if ncold == 0:
                cold = copy.deepcopy(rs)  # initial spectrum is one just read
                starttime = cold.utc
                endtime = starttime
                ncold = 1
            # sums are weighted by durations
                crossZero = False
                crossZeroLat = False
                crossZeroRa = False
                firstlon = rs.gallon # try to keep track of crossing zero in angular coordinates
                firstra = rs.ra #
                cold.ydataA = rs.ydataA * cold.durationSec
            # keep track of observing time for weighted sum
                timesum = rs.durationSec
            else: # else not enough time yet, average cold data
                # fix wrap of longitudes
                cold.count = cold.count + rs.count
                ncold = ncold + 1
                cold.ydataA = cold.ydataA + (rs.ydataA * rs.durationSec)
                cold.ra = cold.ra + (rs.ra * cold.durationSec)
                cold.dec = cold.dec + (rs.dec * cold.durationSec)
                cold.gallon = cold.gallon + (rs.gallon * rs.durationSec)
                cold.gallat = cold.gallat + (rs.gallat * rs.durationSec)
            # keep track of observing time for weighted sum
                endtime = rs.utc
                timesum = timesum + rs.durationSec
            # end if not a enough time
        # end if a cold file
    #end for all files to sum

    if ncrossings > 0:
        crossingsum = crossingsum/float(ncrossings)
    return crossingsum
# check if not data match criteria

# first read through all data and find hot load
names = sys.argv[2:]
names = sorted(names)

minazoffset = -40.
maxazoffset = 40.
mineloffset = -40.
maxeloffset = 40.
eloffsetdelta = maxeloffset-mineloffset
azoffsetdelta = maxazoffset-minazoffset
azoffsetdelta = azoffsetdelta/2.
eloffsetdelta = eloffsetdelta/2.
maxsum = 0.
maxel = 0.
maxaz = 0.
lastmaxsum = 0.
lastmaxel = 0.
lastmaxaz = 0.

gain, vel, lowlatnames = findhotgain( names, avetimesec)

print 'Number of files used for alignment: ',len(lowlatnames)

# iteratively search for maximum
for jjj in range( 7):
    print '%2d: Searching range: %7.2f,%7.2f to %7.2f,%7.2f ' % (jjj,minazoffset,mineloffset, maxazoffset,maxeloffset)

    daz = minazoffset

# for all azimuth offsets
    for naz in range( 5):
        delta = mineloffset
# for all elevation offsets
        for nel in range( 5):

            zerovalue = findzero( lowlatnames, gain, vel, daz, delta, avetimesec)
#        print 'daz, delta, zerovalue:', daz,delta, zerovalue
            sys.stdout.write('\rAt %7.2f, %7.2f: %10.0f \r' % (daz,delta, zerovalue))
            sys.stdout.flush()

# if a better alignment with the galactic plane
            if zerovalue > maxsum:
#                print lastmaxaz, lastmaxel, lastmaxsum
                if lastmaxsum > 0.:
                    maxsum = (lastmaxsum+zerovalue)/2.
                    maxel = ((lastmaxel*lastmaxsum) + (delta*zerovalue))/(lastmaxsum+zerovalue)
                    maxaz = ((lastmaxaz*lastmaxsum) + (daz*zerovalue))/(lastmaxsum+zerovalue)
                else:
                    maxel = delta
                    maxaz = daz
                    maxsum = zerovalue
                lastmaxel = maxel
                lastmaxaz = maxaz
                lastmaxsum = maxsum
                print "   Better Alignment  dAz, dEl: %7.2f,%7.2f; sum: %7.1f" % (lastmaxaz, lastmaxel, maxsum)
            delta = delta + ((maxeloffset-mineloffset)/5.)
        daz = daz + ((maxazoffset-minazoffset)/5.)
            
    print "%2d:  Best Alignment  dAz, dEl: %7.2f,%7.2f; sum: %7.1f" % ( jjj, maxaz, maxel, maxsum)

# reduce range and find max near last max
    azoffsetdelta = azoffsetdelta/2.
    eloffsetdelta = eloffsetdelta/2.

    mineloffset = maxel - eloffsetdelta
    maxeloffset = maxel + eloffsetdelta

    minazoffset = maxaz - azoffsetdelta
    maxazoffset = maxaz + azoffsetdelta
    
