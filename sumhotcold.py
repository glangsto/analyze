#Python Script to compute the integral of the the intensities
#plot the raw data from the observation
#HISTORY
#19Apr22 GIL Update default frequency and velocity ranges
#18Feb16 GIL the Hot and Cold spectra are precommuted for a map
#17Sep22 GIL enable/disable plotting
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
from scipy.signal import savgol_filter
import interpolate
import Moment
import os.path

avetimesec = 3600.
dy = -1.
linelist = [1420.0, 1419.0, 1418.0]  # RFI lines in MHz
linewidth = [11, 5, 9]

#simple argument parsing logic initialization
nargs = len(sys.argv)
timearg = 1
firstfilearg = 2

if nargs < 2:
    print "SUMHOTCOLD: Summarize a set of observations, computing the intensity integrals"
    print "usage: SUMHOTCOLD [plot] <aveTimeSeconds> <hot file> <cold file> <file names>"
    print "where plot           optional string indicating the average spectra should be plotted"
    print "where hot  file      name of the calibration hot  file for this observation"
    print "where cold file      name of the calibratino cold file for this observation"
    print "where aveTimeSeconds time range of data to be averaged before output, in seconds"
    print "where file names     list of files to examine. Should not include hot files"
    print ""
    exit()

# assume no plotting
doPlot = False
firstArg = sys.argv[timearg].upper()
if firstArg == "PLOT":
    doPlot = True
    firstfilearg = firstfilearg + 1
    timearg = timearg + 1

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1e8
xallmax = -9.e9
xallmin = 9.e9
yallmax = -9.e9
yallmin = 9.e9
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
nhot = 0         # number of obs with el < 0
nhigh = 0        # highest galactic latitude
lastfreq = 0.
lastbw = 0.
lastgain = 0.
lastel = 0.
lastaz = 0.
firstaz = -500. # flage no azimuth yet  read
firstdate = ""
lastdate = ""
firstrun = True

#first argument is the averaging time in seconds

avetimesec = float(sys.argv[timearg])
print "Average time: ", avetimesec, " (seconds)"

# hot load file and cold load files are first two file names
hotfilename = sys.argv[firstfilearg]
coldfilename = sys.argv[firstfilearg+1]
# now skip first two files
firstfilearg = firstfilearg+2
names = sys.argv[firstfilearg:]
names = sorted(names)
####################### Set the velocity baseline and averaging range
minvel = -160.
maxvel = 160.
#######################

eloffset = -2.75
azoffset = -34.25
# latest offset after averaging previous best with current best
eloffset = 0.77
azoffset = -34.04
eloffset = 0.0
azoffset = 0.0
# set maximum expected normal System temperature
# this does auto flagging of horn that's tipped over
Tmax = 400.

def fit_baseline( xs, ys, imin, imax):

    xfit = np.concatenate( (xs[imin-20:imin+20],xs[imax-20:imax+20]))
    yfit = np.concatenate( (ys[imin-20:imin+20],ys[imax-20:imax+20]))
# calculate polynomial
    z = np.polyfit(xfit, yfit, 3)
    f = np.poly1d(z)

# calculate new x's and y's
    yout = f(xs)

    return yout

def interpolate_range( minvel, maxvel, vel, ys):
    """
    Interpolate arrays over a range of velocities
    minvel and maxvel are the ranges to interpolate
    input vel array of velocities 
    input ys values to interpolate
    """

    nData = len(vel)
    iref = int(nData/2)
    vref = vel[iref]
    dv   = (vel[iref+2]-vel[iref-2])/4.
    imin = int(((minvel - vref)/dv) + iref)
    imax = int(((maxvel - vref)/dv) + iref) + 1
    if imin > imax:
        temp = imin
        imin = imax
        imax = temp
    if imin < 10:
        print 'Imin Error computing baseline: ', imin
        imin = 10
        if imin >= nData-10:
            print 'Imin Error computing baseline: ', imin
            imin = nData-10
            
    if imax < 11:
        print 'Imax Error computing baseline: ', imax
        imax = 11
    if imax >= nData-11:
        print 'Imax Error computing baseline: ', imax
        imax = nData-11
        
    ya = np.median(ys[(imin-10):(imin+10)])
    yb = np.median(ys[(imax-10):(imax+10)])
    yout = copy.deepcopy(ys)
    slope = (yb-ya)/(imax-imin)
    #    print 'IR: %3d,%3d, ya,b: %7.3f,%7.3f slope: %9.4f' % (imin,imax, ya, yb, slope)
    # finally interpolate over values
    for iii in range( imin, imax):
        yout[iii] = (ya + (slope*(iii-imin)))
    
    return yout

def write_sum( outname, minvel, maxvel, ra, dec, gallon, gallat, time, vel, tsky, el, n):
    """ 
    write out the integral of an array over a relevant velocity range
    """
    global firstrun
    nvel = len(vel)
    nsky = len(tsky)
    if (nsky > nvel):
        print "Error writing sum: Velocites > Intensity array sizes: ", nvel, nsky

    if nvel < 5:
        print 'Not enough values in arrays to compute reliable sums'

    # compute reference velicity and increments
    nvel2 = int(nvel/2)
    dv = (vel[nvel2+2]-vel[nvel2-2])/4.

    vref = vel[nvel2]
    imin = int(((minvel-vref)/dv) + nvel2)
    if imin < 0 or imin > nvel-1:
        print 'Array does not contain minimum velocity: ',minvel
        print 'vel[0]: ', vel[0], 'vel[n-1]: ', vel[nvel-1]
        if imin < 0:
            imin = 0
        else:
            imin = nvel-1

    imax = int(((maxvel-vref)/dv) + nvel2) + 1
    if imax < 0 or imax >= nvel:
        print 'Array does not contain maximum velocity: ',maxvel
        print 'vel[0]: ', vel[0], 'vel[n-1]: ', vel[nvel-1]
        if imax < 0:
            imax = 0
        else:
            imax = nvel-1
    if imax < imin:  # assuming imin > imax for array indexing
        temp = imin
        imin = imax
        imax = temp

    if firstrun:
        print 'Velocity resolution : %9.3f km/sec' % dv 
        print 'Velocity range: %9.2f to %9.2f km/sec' % (minvel, maxvel)
        print 'Index range: %5d to %5d ' % (imin, imax)
        firstrun = False

#    print 'Index range for velocities: ', imin, imax
#    print 'Velocities  for indicies  : ', vel[imin], vel[imax]
#    print 'Temperatures for indicies : ', tsky[imin], tsky[imax]
# first get the RMS in a region not expected to have signal
    vs = vel[imin-20:imin]
    ts = tsky[imin-20:imin]
    (ave,adev,rms,var,skew,curt) = Moment.moment(vs, ts)    

# now get the sum over the range of interest
    vs = vel[imin:imax]
    ts = tsky[imin:imax]
    (ave,adev,sdev,var,skew,curt) = Moment.moment(vs, ts)    
    tmax = max(ts)

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

    
    if not os.path.isfile(outname):
        f = open( outname, 'w')
#        f.write( '#  LON     LAT      T_sum     T_std     V_ave     V_std    Time \n')
#       f.write( '#   RA     DEC     LON     LAT      T_sum     T_std     V_ave     V_std  N  Time \n')
        f.write( '#   RA     DEC     LON     LAT      T_sum +/- T_std   T_max  V_ave+/-V_std  El    N  Time \n')
#        f.write( '#   RA     DEC     LON     LAT      T_sum     T_max     T_std   V_ave  V_std  El    N  Time \n')
#        f.write( '#   RA     DEC      T_sum     T_std     V_ave     V_std    Time \n')
    else:
        f = open( outname, 'a+')

# File = 2017-01-03.sum
# RA  DEC   T_sum  T_std   V_ave  V_std   Time  
# 54.56  -35.32    -5.587     3.039  -624.929   111.851 20:02:38
    tsum = 0
    for iii in range(imax-imin):
        tsum = tsum+ts[iii]
# take absolute value of velocity
    if dv < 0.:
        dv = -dv
    tsum = tsum * dv  # scale integral by velocity resolution
    sdev = adev * dv * np.sqrt(float(imax-imin)) # scale error estimate by approriate number of channels
    print 'Sum, Max:%7.2f %7.2f RMS: %5.2f  %s' % (tsum, tmax, sdev, time)
    if tsum > 0:
        vsdev = rms * dv * np.sqrt(float(imax-imin)*abs(sdev/tsum)) # scale error estimate by approriate number of channels
    label = '%7.2f %7.2f %7.2f %7.2f %9.2f %9.2f %7.2f %7.2f %6.2f %4.1f %3d %s\n' % \
        (ra, dec, gallon, gallat, tsum, sdev, tmax, vave, vsdev, el, n, time)
    f.write(label)
    f.close()

lastgain = 0.

print 'Applying Az,El offset: %6.1f,%6.1f d' % (azoffset, eloffset)
hot = radioastronomy.Spectrum()
print 'Reading hot spectrum: ',hotfilename
hot.read_spec_ast(hotfilename)
hot.telel = hot.telel + eloffset
hot.telaz = hot.telaz + azoffset
#hot.azel2radec()    # compute ra,dec from az,el

nhot = 1

xv = hot.xdata * 1.E-6  # convert frequencies to MHz
yv = hot.ydataA
n = len(yv)
n2 = n/2
yv[n2-1] = ((3.*yv[n2-2])+yv[n2+2])/4.
yv[n2] = (yv[n2-2]+yv[n2+2])/2.
yv[n2+1] = ((yv[n2-2])+(3.*yv[n2+2]))/4.

nData = len(xv)
#previously just copy
#hv = yv
hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

vel = np.zeros(nData)
cold = radioastronomy.Spectrum()
print 'Reading Cold spectrum: ',coldfilename
cold.read_spec_ast(coldfilename)
cold.telel = cold.telel + eloffset
cold.telaz = cold.telaz + azoffset
#cold.azel2radec()    # compute ra,dec from az,el

yv = cold.ydataA
xv = hot.xdata * 1.E-6  # convert frequencies to MHz

cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

# finally compute gain on a channel by chanel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)
    vel[iii] = c * (nuh1 - xv[iii])/nuh1

# now interpolate over galactic velocities to remove galactic HI emission
gain = interpolate_range( minvel, maxvel, vel, gain)

if doPlot:
    fig, ax1 = plt.subplots(figsize=(10, 6))

ymin = 1000.  # init min,max to large values
ymax = 0.
yallmin = ymin
yallmax = ymax

trx = np.zeros(nData)
for iii in range(nData):
    trx[iii] = (cv[iii]/gain[iii]) - tcold


Tsys = np.median(trx)
print "Median Receiver + Antenna Temp: %.1f K" % ( Tsys)

# now recompute gain using only hot load and Trx
# hanning smooth
for iii in range(nData):
    if iii == 0 or iii >= nData-1:
        gain[iii] = hv[iii]/(thot + Tsys)
    else:
        gain[iii] = 0.25*(hv[iii-1]+(2.*hv[iii])+hv[iii+1])/(thot + Tsys)

avetime = datetime.timedelta(seconds=avetimesec)

ncold = 0
nread = 0        
print 'Processing first file: ', names[0]
# now read through all data 
for filename in names:

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    extension = parts[1]

# exclude hot load data for averaging
    if extension == 'hot':
        continue

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    
    dt = str(rs.utc) # get date time
#   separate datetime into date and time parts
    parts = dt.split(' ')
    date = parts[0]
    nd = len(date)
    date = date[2:nd]
    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format
    parts = time.split('.') # trim off fractions of a secon
    time = parts[0]

# update for anly pointing corrections
    rs.telel = rs.telel + eloffset
    rs.telaz = rs.telaz + azoffset
    rs.azel2radec()    # compute ra,dec from az,el

    if firstdate == "":
        firstdate = date

# if first time reading data, set obs parameters
    if lastfreq == 0.:
        lastfreq = rs.centerFreqHz 
        lastbw = rs.bandwidthHz
        lastgain = rs.gains[0]
        firstel = rs.telel
        firstaz = rs.telaz
        cold = copy.deepcopy( rs)
        ncold = 0
        timesum = 0.

    if ncold > 1:
        # time difference is between mid-points of integrations
        dt = rs.utc - cold.utc 
        # add the time since midpoint of latests
        dt = dt + datetime.timedelta(seconds=rs.durationSec/2.)
        # plus time before start of the first
        dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)

        lastdate = date

        newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
        newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl
        
        # if time to average (or end of all files)
        if (dt > avetime) or (filename == sys.argv[nargs-1]) or newObs:
            cold.ydataA = cold.ydataA/float(timesum)
            # have complicated steps to simple get the average time.
            deltatime = endtime - starttime
            # now add offset, in seconds, of half the differencen
            middletime = starttime + datetime.timedelta(seconds=(deltatime.seconds/2.))
            #                print starttime, middletime, middletime
            # now convert to a string and trim out fractions of a second.
            middlestr = str(middletime)
            idot = middlestr.index('.')
            midtime = middlestr[:idot]
    #                print 'Midtime: ', midtime
            # now need to compute average coordinates
            gallon = cold.gallon/float(timesum)
            if gallon > 360.:
                gallon = gallon - 360.
            elif gallon < 0.:
                gallon = gallon + 360.
            ra = cold.ra/float(timesum)
            dec = cold.dec/float(timesum)
            if ra > 360.:
                ra = ra - 360.
            elif ra < 0.:
                ra = ra + 360.
            gallat = cold.gallat/float(timesum)
            az = cold.telaz
            el = cold.telel
# prepare to recalculate galllont
            cold.ra = ra
            cold.dec = dec
            cold.radec2gal()
            gallat = cold.gallat
            gallon = cold.gallon

            xv = cold.xdata * 1.E-6
            yv = cold.ydataA 
##                yv = cold.foldfrequency() * scalefactor
            yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
            xmin = min(xv)
            xmax = max(xv)
            xallmin = min(xmin, xallmin)
            xallmax = max(xmax, xallmax)
            count = cold.count
            note = cold.noteA
                    #print('%s' % note)
            ncolor = min(nmax-1, nplot) 

            tsys = np.zeros(nData)    # initialize arrays
            vel = np.zeros(nData)
            tsky  = np.zeros(nData)    # initialize arrays
            for jjj in range (0, nData):
                vel[jjj] = c * (nuh1 - xv[jjj])/nuh1
                tsys[jjj] = yv[jjj]/gain[jjj]
                tsky[jjj] = tsys[jjj] - Tsys

            Tscan = np.median(tsys)
            if Tscan > Tmax:
                print 'Scan average Temperature is greater than maximum expected: ',Tscan,Tmax
                print 'Scan Time: ',midtime
                ncold = 0
                continue
            label = 'Lon,Lat=%5.1f,%5.1f' % (gallon, gallat)

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
                    
            ymed = np.median(tsky[imin:imax])
            yrms = np.std(tsky[imin-10:imin])
            yrms = yrms + np.std(tsky[imax:imax+10])
            yrms = yrms*.5
#            ya = np.median(tsky[(imin-10):(imin+10)])
#            yb = np.median(tsky[(imax-10):(imax+10)])
#            slope = (yb-ya)/(imax-imin)
#            baseline = tsky
#            print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
# subtract ab baselien
#            for iii in range( nData):
#                baseline = (ya + (slope*(iii-imin)))
            baseline = interpolate_range( minvel, maxvel, vel, tsky)
#            baseline = fit_baseline( vel, tsky, imin-20, imax+20)
            tsky = tsky - baseline

# no longer folding so no factor of two needed
#                tsky = 2. * tsky

            ymin = min(tsky[(2*nData/4):(3*nData/4)])
            ymax = max(tsky[(2*nData/4):(3*nData/4)])
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            label = '%s Az,El: %5s,%5s, Lon,Lat: %5.1f,%5.1f' % (time, az, el, gallon, gallat)
#                label = '%s Lon,Lat=%5.1f,%5.1f' % (colddate+' '+time,gallon,gallat)
            label = '%s Lon,Lat=%5.1f,%5.1f' % (midtime, gallon, gallat)
#                timesum = cold.utc + rs.utc
#                timesum = timesum/2.
            hours = time[0:2]
            outname = date + '_' + hours + '.sum'
            # compute sums over a smaller range that the fit
            write_sum( outname, minvel+10., maxvel-10., ra, dec, gallon, gallat, date + ' ' + time, vel, tsky, el, ncold)

            if ncold > 1:
                print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, yrms, ncold, label)
                if doPlot:
                    if gallat < 7.5 and gallat > -7.5:
                        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
                    elif gallat < 15. and gallat > -15.:
                        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
                    else:
                        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
                    nplot = nplot + 1

                ncold = 0

        # end if a new observation
        if newObs:
            print "Change in observing parameters: "
            if lastfreq != rs.centerFreqHz:
                print "LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz"
                lastfreq = rs.centerFreqHz
            if lastbw != rs.bandwidthHz:
                print "LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz"
                lastbw = rs.bandwidthHz
            if lastgain != rs.gains[0]:
                if lastgain != 0:                    
                    print "LastGain: ", lastgain, "New: ", rs.gains[0], " dB"
                lastgain = rs.gains[0]
            if newAzEl:
                if lastaz != 0 or lastel != 0:
                    print "LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees"
                lastaz = rs.telaz
                lastel = rs.telel

# if this was a new obs; restart the sums
    if ncold == 0:
        cold = copy.deepcopy(rs)  # initial spectrum is one just read
        colddate = date
        starttime = cold.utc
        endtime = starttime
        ncold = 1
        #            print 'Xmin: ', min(cold.xdata)/1e6, 'Xmax: ', max(cold.xdata),' MHz'
        # sums are weighted by durations
        crossZero = False
        crossZeroRa = False
        firstlon = rs.gallon # try to keep track of crossing zero in angular coordinates
        firstra = rs.ra #
        cold.ydataA = rs.ydataA * cold.durationSec
        cold.gallat = rs.gallat * cold.durationSec
        cold.gallon = rs.gallon * cold.durationSec
        cold.ra = rs.ra * cold.durationSec
        cold.dec = rs.dec * cold.durationSec
        # keep track of observing time for weighted sum
        timesum = rs.durationSec
    else: # else not enough time yet, average cold data
        # fix wrap of longitudes
        if abs(rs.gallon - firstlon) > 180:
            crossZero = True
            if rs.gallon > firstlon:
                rs.gallon = rs.gallon - 360.
            else:
                rs.gallon = rs.gallon + 360.
            if abs(rs.ra - firstra) > 180:
                #                print 'Zero RA Crossing: ', rs.ra, firstra
                crossZeroRa = True
            if rs.ra > firstra:
                rs.ra = rs.ra - 360.
            else:
                rs.ra = rs.ra + 360.
            if abs(rs.ra - firstra) > 180:
                if crossZeroRa:
                    print "A problem with this ra", rs.ra," compared with ",firstra
                    continue # jump out without updateing sums
        cold.count = cold.count + rs.count
        ncold = ncold + 1
        cold.ydataA = cold.ydataA + (rs.ydataA * cold.durationSec)
        cold.ra = cold.ra + (rs.ra * cold.durationSec)
        cold.dec = cold.dec + (rs.dec * cold.durationSec)
        cold.gallon = cold.gallon + (rs.gallon * cold.durationSec)
        cold.gallat = cold.gallat + (rs.gallat * cold.durationSec)
    # keep track of observing time for weighted sum
        endtime = rs.utc
        timesum = timesum + rs.durationSec
        # end if not a enough time
    # end if a cold file
#end for all files to sum

# check if not data match criteria

if nplot < 1:
    exit()

if not doPlot:
    print "Computed ",nplot," Integrated Intensities"
    exit()

if firstdate != lastdate:
    date = firstdate + " - " + lastdate
else:
    date = firstdate

if firstaz < -360.:
    firstaz = rs.telaz
    firstel = rs.telel

mytitle = "%s    Az=%6.1f, El=%6.1f" % (date, firstaz, firstel)
fig.canvas.set_window_title(mytitle)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
#plt.xlim(-250., 300.)
#plt.xlim(-300., 600.)
plt.xlim( 2.*minvel, 2*maxvel)
#plt.xlim(xallmin,xallmax)
#plt.xlim(-120., 130.)
plt.ylim(yallmin*.97, 1.15*yallmax)
#plt.ylim(-5., 80.)
plt.title(mytitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
    
