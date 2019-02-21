#Python Script to plot average NSF spectra.
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
from PIL import Image
from pylab import *
import sys
import datetime
import statistics
import radioastronomy
import copy

def raDecToPix( ra, dec, nx, ny, transform):
    """
    Transform astronomical coordinate to x,y pixels of image
    Inputs:
    ra = right ascension in hours
    dec = declination in degrees
    nx = image x size in pixels
    ny = image y size in pixels
    transform = coordinate transform parameters (dummy)
    """

    if transform == 1:
# set the corners of the image for placement of spectra locations
        x0 = 360.  # hours
        xz = 0.   # hours
        y0 = 90.  # degrees
        yz = -90. # degrees

        dx = xz - x0
        dy = yz - y0
        dra = ra - x0
        ddec = dec - y0
#        print 'dx, dy, dra, ddec: ', dx, dy, dra, ddec
    
        ix = int(nx * dra / dx)
        iy = int(ny * ddec / dy)

        print 'ix, iy: ', ix, iy
    else:
        ix = int(nx/2)
        iy = int(ny/2)
# end of raDecToPix
    return ix, iy
 
avetimesec = 120.
dy = -1.

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
xa = 150
xb = 1024-xa

#for symbol, value in locals().items():
#    print symbol, value

nplot = 0
nhot = 0
ncold = 0
lastfreq = 0.
lastbw = 0.
lastgain = 0.
lastel = 0.
lastaz = 0.
firstdate = ""
lastdate = ""

#first argument is the averaging time in seconds

avetimesec = float(sys.argv[1])
print "Average time: ", avetimesec, " (seconds)"

# first read through all data and find hot load
names = sys.argv[2:]
names = sorted(names)
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el

    if rs.telel < 0:
        if nhot == 0:
            hot = copy.deepcopy( rs)
            nhot = 1
        else:
            hot.ydataA = hot.ydataA + rs.ydataA
            hot.count = hot.count + rs.count
            nhot = nhot + 1
    else:
        ncold = ncold + 1

if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA / float(nhot)
    print "Found %3d hot load obs" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.E-6
yv = hot.ydataA
nData = len(xv)
hv = yv

vel = np.zeros(nData)
for jjj in range (0, nData):
    vel[jjj] = c * (xv[jjj] - nuh1)/nuh1

# will need smoothed hot load values in remaining calc
#for iii in range(1,(nData-2)):
#    hv[iii] = (yv[iii-2]+yv[iii-1]+yv[iii]+yv[iii+1]+yv[iii+2])/5.
hv = yv

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
# make arrays to hold coordinates
nobs = ncold
print 'NCold sky obs: ',ncold
ras = np.zeros(100000)
decs = np.zeros(100000)
obs = np.zeros(100000)
vmax = np.zeros(100000)
print ras[1], decs[1], obs[1]
print ras[0], decs[0], obs[0]

ncold = 0
#print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
#plt.plot(xv, yv, colors[0], linestyle=linestyles[0],label=label)

# condition data for avoid divide by zero later (2 depends on scalefactor)
for iii in range(nData):
    if hv[iii] < 2.:
        hv[iii] = 2.

avetime = datetime.timedelta(seconds=avetimesec)

# map in galactic coordinates
filename = "/Users/glangsto/Desktop/Research/milkyway_21cm.jpg"
# map in ra-dec coordinates
filename = "/Users/glangsto/Desktop/Research/HEC_visible_sky_location_HR_crop.png"
filename = "/Users/glangsto/Desktop/Research/HEC_visible_sky_location_crop.png"
filename = "/Users/glangsto/Desktop/Research/1000Mess.png"

# read image to array
plt.tight_layout()
im = array(Image.open(filename))
nx = im.shape[1]
ny = im.shape[0]
#print im.shape, im.dtype
xticks = [0, nx/4, nx/2, 3*nx/4, nx-1]
yticks = [0, ny/4, ny/2, 3*ny/4, ny-1]


ax1.set_yticks(yticks)
ax1.set_xticks(xticks)
xlabels = ax1.set_xticklabels(('24h','18h','12h','6h','0h'))
ylabels = ax1.set_yticklabels(('90d','45d','0d','-45d','-90d'))
ax1.set_xlabel("Local Siderial Time (hours)")
ax1.set_ylabel("Declination (degrees)")

#, extent=[24,0,-90,90])
imshow(im)


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
    if firstdate == "":
        firstdate = date
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

# if first time reading data, set obs parameters
        if lastfreq == 0.:
            lastfreq = rs.centerFreqHz 
            lastbw = rs.bandwidthHz
            lastgain = rs.gains[0]
            cold = copy.deepcopy( rs)
            crosszero = False
            firstlon = rs.gallon
            firstra = rs.ra
            ncold = 0

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
            if newObs:
                print "Change in observing parameters: "
                if lastfreq != rs.centerFreqHz:
                    print "LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz"
                    lastfreq = rs.centerFreqHz
                if lastbw != rs.bandwidthHz:
                    print "LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz"
                    lastbw = rs.bandwidthHz
                if lastgain != rs.gains[0]:
                    print "LastGain: ", lastgain, "New: ", rs.gains[0], " dB"
                    lastgain = rs.gains[0]
                if newAzEl:
                    print "LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees"
                    lastaz = rs.telaz
                    lastel = rs.telel

           # if time to average (or end of all files)
            if (dt > avetime) or (filename == sys.argv[nargs-1]) or newObs:
                cold.ydataA = cold.ydataA/float(timesum)

                gallon = cold.gallon/float(timesum)
                if gallon > 360.:
                    gallon = gallon - 360.
                elif gallon < 0.:
                    gallon = gallon + 360.

                ra = cold.ra/float(timesum)
                if ra > 360.:
                    ra = ra - 360.
                elif ra < 0.:
                    ra = ra + 360.
                dec = cold.dec/float(timesum)

#                print nplot, ra, dec
                ras[nplot] = ra
                decs[nplot] = dec
                    
                gallat = cold.gallat/float(timesum)
                az = cold.telaz
                el = cold.telel

                xv = cold.xdata * 1.E-6
                yv = cold.ydataA * scalefactor

                xmin = min(xv)
                xmax = max(xv)
                xallmin = min(xmin, xallmin)
                xallmax = max(xmax, xallmax)
                count = cold.count
                note = cold.noteA
                    #print('%s' % note)
                ncolor = min(nmax-1, nplot) 

                tsys = np.zeros(nData)    # initialize arrays
                iss = np.zeros(nData)
                Z = np.zeros(nData)
                oneMZ = np.zeros(nData)
                for jjj in range(0, nData):
                    iss[jjj] = jjj
                    if yv[jjj] < 1.:
                        yv[jjj] = 1.

                    Z[jjj] = yv[jjj]/hv[jjj]
                    oneMZ[jjj] = 1. - Z[jjj]
                    if oneMZ[jjj] < .001:
                        oneMZ[jjj] = 0.001
                    tsys[jjj] = ((Z[jjj]*thot) - tcold)/oneMZ[jjj]

                vel = np.zeros(nData)
                for jjj in range (0, nData):
                    vel[jjj] = c * (xv[jjj] - nuh1)/nuh1

                tsky  = np.zeros(nData)    # initialize arrays
                S     = np.zeros(nData)    # initialize arrays

                for jjj in range (0, nData):
                    vel[jjj] = c * (xv[jjj] - nuh1)/nuh1
                    S[jjj] = hv[jjj]/(tsys[jjj] + thot)
                    tsky[jjj] = yv[jjj]/S[jjj]

                # remove spike in center of the plot
                icenter = 512
                tsky[icenter] = (tsky[icenter-2] + tsky[icenter+2])*.5
                tsky[icenter-1] = (3.*tsky[icenter-2] + tsky[icenter+2])*.25
                tsky[icenter+1] = (tsky[icenter-2] + 3.*tsky[icenter+2])*.25

                label = 'Lon,Lat=%5.1f,%5.1f' % (gallon, gallat)

                ymed = statistics.median(tsky)
                ya = statistics.median(tsky[(xa-10):(xa+10)])
                yb = statistics.median(tsky[(xb-10):(xb+10)])
                slope = (yb-ya)/(xb-xa)
#                baseline = tsky
#                print 'ya,b: %6.1f,%6.1f; slope: %8.4f' % (ya, yb, slope)
                imax = 0
                for iii in range( 0, nData):
                    tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))
                    if tsky[iii] > tsky[tmax]:
                        imax = iii
                obs[nplot] = tsky[imax]
                vmax[nplot] = vel[imax]

                ymin = min(tsky[(nData/8):(7*nData/8)])
                ymax = max(tsky[(nData/8):(7*nData/8)])
                yallmin = min(ymin,yallmin)
                yallmax = max(ymax,yallmax)
                label = '%s, Az,El: %5s,%5s, Lon,Lat: %5.1f,%5.1f' % (time, az, el, gallon, gallat)
                label = '%s Lon,Lat=%5.1f,%5.1f' % (time,gallon,gallat)
                print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)
#                if ncold > 2:
#                    plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
#                plt.plot(vel, yv, colors[ncolor], linestyle=linestyles[ncolor],label=label)

                nplot = nplot + 1
                ncold = 0
            # end if a new observation
# if this was a new obs; restart the sums
        if ncold == 0:
            cold = copy.deepcopy(rs)  # initial spectrum is one just read
            ncold = 1
#            print 'Xmin: ', min(cold.xdata)/1e6, 'Xmax: ', max(cold.xdata),' MHz'
            # sums are weighted by durations
            firstlon = rs.gallon
            firstra  = rs.ra
            cold.ydataA = rs.ydataA * rs.durationSec
            cold.gallat = rs.gallat * rs.durationSec
            cold.gallon = rs.gallon * rs.durationSec
            cold.ra = rs.ra * rs.durationSec
            cold.dec = rs.dec * rs.durationSec
            # keep track of observing time for weighted sum
            timesum = cold.durationSec
        else: # else ont enough time yet, average cold data
            cold.count = cold.count + rs.count
            ncold = ncold + 1
            cold.ydataA = cold.ydataA + (rs.ydataA * rs.durationSec)
            # fix wrap of longitudes
            if abs(rs.gallon - firstlon) > 180:
                if rs.gallon > firstlon:
                    rs.gallon = rs.gallon - 360.
                else:
                    rs.gallon = rs.gallon + 360.
            if abs(rs.ra - firstra) > 180:
                if rs.ra > firstra:
                    rs.ra = rs.ra - 360.
                else:
                    rs.ra = rs.ra + 360.
            cold.ra = cold.ra + (rs.ra * rs.durationSec)
            cold.dec = cold.dec + (rs.dec * rs.durationSec)
            cold.gallon = cold.gallon + (rs.gallon * rs.durationSec)
            cold.gallat = cold.gallat + (rs.gallat * rs.durationSec)
            # keep track of observing time for weighted sum
            timesum = timesum + cold.durationSec
#                plt.plot(iss, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            # end if not a enough time
        # end if a cold file
    #end for all files to sum

#plt.xlim(xallmin,xallmax)
if firstdate != lastdate:
    date = firstdate + " - " + lastdate
else:
    date = firstdate

# some points

ras = ras[0:nplot]
decs = decs[0:nplot]
obs = obs[0:nplot]
vmax = vmax[0:nplot]

ixs = ras
iys = decs
transform = 1

# for all average coordinates
for iii in range(nplot):
    ixs[iii],iys[iii] = raDecToPix( ras[iii], decs[iii], nx, ny, transform)

# plot the points with red star-markers
plot(ixs, iys,'r*')

# add title and show the plot
title('Plotting: %s' % (filename))
show()

mytitle = "%s    Az=%6.1f, El=%6.1f" % (date, az, el)
#fig.canvas.set_window_title(mytitle)
#for tick in ax1.xaxis.get_major_ticks():
#    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
#for tick in ax1.yaxis.get_major_ticks():
#    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')
#plt.xlim(-250., 300.)
#plt.ylim(yallmin*.97, 1.15*yallmax)
#plt.ylim( 0, max(hv))
#plt.title(mytitle, fontsize=16)
#plt.xlabel('Velocity (km/sec)', fontsize=16)
#plt.ylabel('Intensity (Kelvins)', fontsize=16)
#plt.legend(loc='upper left')
#plt.legend(loc='upper right')
#plt.show()
