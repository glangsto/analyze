#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
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
import statistics
import radioastronomy
import copy
from scipy.signal import savgol_filter
import interpolate
import Moment
import os.path

avetimesec = 3600.
dy = -1.
linelist = [1420.0, 1419.0, 1418.0]  # RFI lines in MHz
linewidth = [9, 5, 9]

nargs = len(sys.argv)

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
ncold = 0
nhigh = 0        # highest galactic latitude
minGlat = +90.
maxGlat = -90.
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

avetimesec = float(sys.argv[1])
print "Average time: ", avetimesec, " (seconds)"

# first read through all data and find hot load
names = sys.argv[2:]
names = sorted(names)
minvel = -160.
maxvel = 160.
minvel = -130.
maxvel = 130.

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
        
    ya = statistics.median(ys[(imin-10):(imin+10)])
    yb = statistics.median(ys[(imax-10):(imax+10)])
    slope = (yb-ya)/(imax-imin)
    # finally interpolate over values
    for iii in range( imin, imax):
        ys[iii] = (ya + (slope*(iii-imin)))
    
    return ys

def write_sum( outname, minvel, maxvel, ra, dec, gallon, gallat, time, vel, tsky, n):
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
    if firstrun:
        print 'Velocity resolution : %9.3f km/sec' % dv 
        print 'Velocity range: %9.2f to %9.2f km/sec' % (minvel, maxvel)
        print 'Index range: %5d to %5d ' % (imin, imax)
        print 'ts: ', ts
        firstrun = False

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
        f.write( '#   RA     DEC     LON     LAT      T_sum     T_std     V_ave     V_std  N  Time \n')
#        f.write( '#   RA     DEC      T_sum     T_std     V_ave     V_std    Time \n')
    else:
        f = open( outname, 'a+')

# File = 2017-01-03.sum
# RA  DEC   T_sum  T_std   V_ave  V_std   Time  
# 54.56  -35.32    -5.587     3.039  -624.929   111.851 20:02:38
    tsum = 0
    for iii in range(imax-imin):
        tsum = tsum+ts[iii]
    tsum = tsum * abs(dv)  # scale integral by velocity resolution
    sdev = rms * abs(dv) * np.sqrt(float(imax-imin)) # scale error estimate by approriate number of channels
    if tsum > 0:
        vsdev = rms * dv * np.sqrt(float(imax-imin)*abs(sdev/tsum)) # scale error estimate by approriate number of channels
    label = '%7.2f %7.2f %7.2f %7.2f %9.2f %9.3f %9.3f %9.3f %3d %s\n' % (ra, dec, gallon, gallat, tsum, sdev, vave, vsdev, n, time)
    f.write(label)
    f.close()

lastgain = 0.
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el

    if rs.telel < 0:
        if rs.gains[0] != lastgain:
            if lastgain != 0:
                print 'Gain Change: ', lastgain, ' -> ', rs.gains[0]
            lastgain = rs.gains[0]
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
yv = hot.foldfrequency()
#yv[511] = ((3.*yv[510])+yv[514])/4.
#yv[512] = (yv[510]+yv[514])/2.
#yv[513] = ((yv[510])+(3.*yv[514]))/4.
nData = len(xv)
#previously just copy
#hv = yv
hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

vel = np.zeros(nData)
for jjj in range (0, nData):
    vel[jjj] = c * (nuh1 - xv[jjj])/nuh1

# will need smoothed hot load values in remaining calc
#for iii in range(1,(nData-2)):
#    hv[iii] = (yv[iii-2]+yv[iii-1]+yv[iii]+yv[iii+1]+yv[iii+2])/5.
#hv = yv

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

# but if no low latitude data are available use minimum
# now average coldest data for calibration
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
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

if nhigh < 1.:
    print "No high galactic latitude data: can not calibrate"
    exit()
else:
    high.ydataA = scalefactor * high.ydataA/nhigh
    print "Found %d High Galactic Latidue spectra" % (nhigh)
    yv = high.ydataA
    yv = high.foldfrequency()
#    yv[511] = ((3.*yv[510])+yv[514])/4.
#    yv[512] = (yv[510]+yv[514])/2.
#    yv[513] = ((yv[510])+(3.*yv[514]))/4.

    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

# finally compute gain on a channel by chanel basis
gain = np.zeros(nData)
for iii in range(nData):
    gain[iii] = (hv[iii] - cv[iii])/(thot - tcold)
    vel[iii] = c * (nuh1 - xv[iii])/nuh1

# now interpolate over galactic velocities
gain = interpolate_range( minvel, maxvel, vel, gain)

fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
ymed = statistics.median(yv)
count = hot.count
ncold = 0

trx = np.zeros(nData)
for iii in range(nData):
    trx[iii] = (cv[iii]/gain[iii]) - tcold

Tsys = statistics.median(trx)
print "Median Receiver + Antenna Temp: ", Tsys

#plt.plot(xv, trx, colors[1], linestyle=linestyles[0],label="Tsys")
#plt.show()

avetime = datetime.timedelta(seconds=avetimesec)

#plt.plot(vel, hv, colors[0], linestyle=linestyles[3],label="hot")

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
    rs.azel2radec()    # compute ra,dec from az,el

# if a sky observation
    if rs.telel > 0.:

#        if rs.gallat < 70.0 and rs.gallat > -70.:
#            continue

        if firstdate == "":
            firstdate = date
        # do not put in the header the elevation of the off source
        if rs.gallat < 40. and rs.gallat > -40.:
            firstel = rs.telel
            firstaz = rs.telaz

# if first time reading data, set obs parameters
        if lastfreq == 0.:
            lastfreq = rs.centerFreqHz 
            lastbw = rs.bandwidthHz
            lastgain = rs.gains[0]
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

                xv = cold.xdata * 1.E-6
                yv = cold.ydataA * scalefactor
                yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
# fold to clear up band flip problem
                cold.ydataA = yv
                yv = cold.foldfrequency()
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
                    
                ymed = statistics.median(tsky[imin:imax])
                ya = statistics.median(tsky[(imin-10):(imin+10)])
                yb = statistics.median(tsky[(imax-10):(imax+10)])
                slope = (yb-ya)/(imax-imin)
#                baseline = tsky
#                print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
# subtract ab baselien
#                for iii in range( 0, nData):
#                    tsky[iii] = tsky[iii] - (ya + (slope*(iii-imin)))
                baseline = fit_baseline( vel, tsky, imin-20, imax+20)
                tsky = tsky - baseline

# folding so factor of two needed
                tsky = 2. * tsky

                ymin = min(tsky[(2*nData/4):(3*nData/4)])
                ymax = max(tsky[(2*nData/4):(3*nData/4)])
                yallmin = min(ymin,yallmin)
                yallmax = max(ymax,yallmax)
                label = '%s, Az,El: %5s,%5s, Lon,Lat: %5.1f,%5.1f' % (time, az, el, gallon, gallat)
#                label = '%s Lon,Lat=%5.1f,%5.1f' % (colddate+' '+time,gallon,gallat)
                label = '%s Lon,Lat=%5.1f,%5.1f' % (midtime, gallon, gallat)
#                timesum = cold.utc + rs.utc
#                timesum = timesum/2.
                hours = time[0:2]
                outname = date + '_' + hours + '.sum'
                write_sum( outname, -90., +90., ra, dec, gallon, gallat, date + ' ' + time, vel, tsky, ncold)

                if ncold > 1:
                    print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)
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
#                plt.plot(iss, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            # end if not a enough time
        # end if a cold file
    #end for all files to sum

# check if not data match criteria

if nplot < 1:
    print 'No Data passed selection critera, no plots'
    exit()

#plt.xlim(xallmin,xallmax)
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
#plt.xlim(-150., 350.)
plt.xlim(-250., 450.)
#plt.xlim(xallmin,xallmax)
#plt.xlim(-120., 130.)
plt.ylim(yallmin*.97, 1.15*yallmax)
#plt.ylim(yallmin*.97, 70)
#plt.ylim( 0, max(hv))
plt.title(mytitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.show()
