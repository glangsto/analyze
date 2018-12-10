#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#18FEB06 GIL only output if nfiles > 1
#17AUG17 GIL Note elevation range
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

EPSILON = 0.0001
avetimesec = 3600.
avetimesec = 900.
#avetimesec = 60.
avetimesec = 120.
dy = -1.

# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False
flagCenter = True
# put list of RFI features here, for interpolation later
flagRfi = False
flagRfi = True
linelist = [1400.00, 1420.00]  # RFI lines in MHz
linewidth = [5, 5] # line width is in channels

nargs = len(sys.argv)

# these are the line colors and styles for plotting
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
minel = 200.
maxel = -200.

#first argument is the averaging time in seconds
timearg = 1
namearg = 2

# if folding data;  This was for problem diagnosis.  No longer used much
if sys.argv[1] == '-f':
    print 'Folding specectra'
    dofold = True
    timearg = timearg+1
    namearg = namearg+1
else:
    dofold = False

avetimesec = float(sys.argv[timearg])
print "Average time: ", avetimesec, " (seconds)"

# first read through all data and find hot load
names = sys.argv[namearg:]
names = sorted(names)

rs = radioastronomy.Spectrum()
for filename in names:

    parts = filename.split('/')
    nparts = len(parts)
    if nparts == 1:
        aname = parts[0]
    else:
        aname = parts[nparts-1]
    parts = aname.split('.')
    nparts = len(parts)
    if nparts < 2:
        print 'File is not an astronomy file: ',filename
        continue
    else:
        extension = parts[nparts-1]
    extension = extension.upper()
    if (extension != 'HOT') and (extension != 'AST') and (extension != 'CLD'):
        print 'Extension not recognized : ', parts[nparts-1]
        continue

    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if dofold:
        rs.foldfrequency()

    if rs.telel < 0:
        if nhot == 0:
            hot = copy.deepcopy( rs)
            nhot = 1
        else:
            hot.ydataA = hot.ydataA + rs.ydataA
            hot.count = hot.count + rs.count
            nhot = nhot + 1
    else:
        if minel > rs.telel:
            minel = rs.telel
        if maxel < rs.telel:
            maxel = rs.telel

if nhot > 0:
    print "Found %3d hot load observations" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.E-6
nData = len( xv)
n6 = int(nData/6)
n56 = 5*n6
yv = hot.ydataA
#yv = hot.foldfrequency()

hotmedian = statistics.median(yv[n6:n56])
if hotmedian > .001:
    scalefactor = 1.0

hot.ydataA = scalefactor * hot.ydataA / float(nhot)

# if flagging RFI 
if flagRfi:
    hv = radioastronomy.lines( linelist, linewidth, xv, yv) # interpolate rfi
else:
    hv = yv

fig, ax1 = plt.subplots(figsize=(10, 6))
#plt.hold(True)
az = hot.telaz
el = hot.telel
ymin = 1000.  # initi to large values
ymax = 0.
yallmin = ymin
yallmax = ymax
ymed = statistics.median(yv)
count = hot.count
ncold = 0

def compute_tsky( xv, yv, hv, thot, tcold):
    if flagRfi:
        yv = radioastronomy.lines( linelist, linewidth, xv, yv) # interpolate rfi
    nData = len(yv)    
    tsys = np.zeros(nData)    # initialize arrays

    Z = np.zeros(nData)
    oneMZ = np.zeros(nData)
    # For full Temp calibration, a spectrum taken at high elevation away from 
    # The galactic plan is used.   For this program the cold spectrum must be
    # the spectrum being calibrated.   See the M command for comparision
    epsilons = np.full( nData, EPSILON)
    yv = np.maximum( yv, epsilons)
    hv = np.maximum( hv, epsilons)
        # comput the cold/hot ratio
    Z = yv/hv
    oneMZ = np.full( nData, 1.) - Z
    oneMZ = np.maximum( oneMZ, epsilons)

    # the cold, receiver, temperature is this function
    tsys = ((Z*thot) - tcold)/oneMZ
    
    n6 = int(nData/6)
    n56 = 5*n6

    tsysmedian = statistics.median( tsys[n6:n56])

    tsky  = np.zeros(nData)    # initialize arrays
    S     = np.zeros(nData)    # initialize arrays

    # The system gain S is computed assuming a tsys is the cold load
    S = np.full( nData, tsysmedian+thot)/hv
    # scale the observed instensity in counts to Kelvins.
    tsky = S*yv

    if flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
        icenter = int(nData/2)
        tsky[icenter] = (tsky[icenter-2] + tsky[icenter+2])*.5
        tsky[icenter-1] = (3.*tsky[icenter-2] + tsky[icenter+2])*.25
        tsky[icenter+1] = (tsky[icenter-2] + 3.*tsky[icenter+2])*.25

    return tsky

# condition data for avoid divide by zero later (2 depends on scalefactor)
hv = np.maximum( hv, np.full( nData, 2.))

avetime = datetime.timedelta(seconds=avetimesec)
dt = datetime.timedelta(seconds=0.)
newObs = False

nRead = 0        
nFiles = len(names)
print nFiles

# now read through all data and average cold sky obs
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if dofold:
        rs.foldfrequency()

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    dtstr = str(rs.utc)
    parts = dtstr.split(' ')
    date = parts[0]
    nd = len(date)
    date = date[2:nd]
    if firstdate == "":
        firstdate = date
    nRead = nRead + 1

# if not a sky observation; skip further processing
    if rs.telel < 0.:
        continue

    time = parts[1]
    time = time.replace('_', ':')  # put time back in normal hh:mm:ss format
    parts = time.split('.')        # trim off fractional seconds
    time = parts[0]

# if first time reading data, set obs parameters
    if lastfreq == 0.:
        lastfreq = rs.centerFreqHz 
        lastbw = rs.bandwidthHz
        lastgain = rs.gains[0]
        lastaz = rs.telaz
        lastel = rs.telel
        cold = copy.deepcopy( rs)
        firstutc = rs.utc
        lastutc = rs.utc
        ncold = 0

    if ncold > 0:
        # time difference is between mid-points of integrations
        dt = rs.utc - cold.utc 
        # add the time since midpoint of latests
        dt = dt + datetime.timedelta(seconds=rs.durationSec/2.)
        # plus time before start of the first
        dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)
        lastutc = rs.utc

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

    if nRead == nFiles:   # if this is the last file, must force output
        newObs = True
                    
        # if time to average (or end of all files)
    if (dt > avetime) or (filename == sys.argv[nargs-1]) or newObs:
        cold.ydataA = cold.ydataA/float(timesum)
            
        gallon = cold.gallon/float(timesum)
        gallat = cold.gallat/float(timesum)
        az = cold.telaz
        el = cold.telel

        # convert to MHz from Hz
        xv = cold.xdata * 1.E-6

        xmin = min(xv)
        xmax = max(xv)
        xallmin = min(xmin, xallmin)
        xallmax = max(xmax, xallmax)
        count = cold.count
        note = cold.noteA
        ncolor = min(nmax-1, nplot) 

        yv = cold.ydataA * scalefactor
        tsky = compute_tsky( xv, yv, hv, thot, tcold)

        ymin = min(tsky[(nData/8):(7*nData/8)])
        ymax = max(tsky[(nData/8):(7*nData/8)])
        yallmin = min(ymin,yallmin)
        yallmax = max(ymax,yallmax)
        ymed = statistics.median(tsky)
        aveutc,duration = radioastronomy.aveutcs( firstutc, lastutc)
        strnow = aveutc.isoformat()
        datestr = strnow.split('.')
        daypart = datestr[0]
        time = daypart[12:19]
        
        if firstdate == lastdate: 
            label = '%s Lon,Lat=%5.1f,%5.1f' % (time, gallon, gallat)
        else:
            label = '%s %s Lon,Lat=%5.1f,%5.1f' % (date, time, gallon, gallat)
            
        if minel == maxel: 
            ellabel = ''
        else:
            ellabel = ' El=%4.1f' % (el)
        label = label + ellabel
        print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)

        # plot thicker lines when near the galactic plane
        if gallat < 7.5 and gallat > -7.5:
            plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
        elif gallat < 15. and gallat > -15.:
            plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
        else:
            plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            
        nplot = nplot + 1
        ncold = 0

# if this was a new obs; restart the sums
    if ncold == 0:
        cold = rs  # initial spectrum is one just read
        ncold = 1
        # sums are weighted by durations
        firstlon = rs.gallon
        cold.ydataA = cold.ydataA * cold.durationSec
        cold.gallat = cold.gallat * cold.durationSec
        cold.gallon = cold.gallon * cold.durationSec
    # keep track of observing time for weighted sum
        timesum = cold.durationSec
    else: # else ont enough time yet, average cold data
        cold.count = cold.count + rs.count
        ncold = ncold + 1
        cold.ydataA = cold.ydataA + (rs.ydataA * cold.durationSec)
        # fix wrap of longitudes
        if abs(rs.gallon - firstlon) > 180:
            crosZero = True
            if rs.gallon > firstlon:
                rs.gallon = rs.gallon - 360.
            else:
                rs.gallon = rs.gallon + 360.
        cold.gallon = cold.gallon + (rs.gallon * cold.durationSec)
        cold.gallat = cold.gallat + (rs.gallat * cold.durationSec)
        # keep track of observing time for weighted sum
        timesum = timesum + cold.durationSec
        # end if not a enough time
        # end if a cold file

# if here, all files read, plot the remaining observation
if ncold > 0:
    cold.ydataA = cold.ydataA/float(timesum)

    gallon = cold.gallon/float(timesum)
    gallat = cold.gallat/float(timesum)
    az = cold.telaz
    el = cold.telel

    # convert to MHz from Hz
    xv = cold.xdata * 1.E-6
    
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin, xallmin)
    xallmax = max(xmax, xallmax)
    count = cold.count
    note = cold.noteA
                    #print('%s' % note)
    ncolor = min(nmax-1, nplot) 

    yv = cold.ydataA * scalefactor
    tsky = compute_tsky( xv, yv, hv, thot, tcold)
    
    ymin = min(tsky[(nData/8):(7*nData/8)])
    ymax = max(tsky[(nData/8):(7*nData/8)])
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    ymed = statistics.median(tsky)
    if firstdate == lastdate: 
        label = '%s Lon,Lat=%5.1f,%5.1f' % (time, gallon, gallat)
    else:
        label = '%s %s Lon,Lat=%5.1f,%5.1f' % (date, time, gallon, gallat)

    if minel == maxel: 
        ellabel = ''
    else:
        ellabel = ' El=%4.1f' % (el)
    label = label + ellabel
    print ' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, ncold, label)

    # plot thicker lines when near the galactic plane
    if gallat < 7.5 and gallat > -7.5:
        plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
    elif gallat < 15. and gallat > -15.:
        plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
    else:
        plt.plot(xv, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)

    nplot = nplot + 1
    ncold = 0

#plt.xlim(xallmin,xallmax)
lastdate = date
# if observations span several days
if lastdate != firstdate:
    date = firstdate + " to " + lastdate

# if change in elevation 
if minel == maxel:
    mytitle = "%s    Az=%6.1f, El=%6.1f" % (date, az, minel)
else:
    mytitle = "%s    Az=%6.1f, El=%6.1f to %6.1f" % (date, az, minel, maxel)

fig.canvas.set_window_title(mytitle)

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
#plt.xlim(-300., 500.)
plt.xlim(xallmin, xallmax)
#plt.xlim(1422.1, 1422.9)
#plt.ylim((yallmin*.97)-1., 1.25*yallmax)
#plt.ylim(yallmin*.97, 400.)
plt.title(mytitle, fontsize=16)
plt.xlabel('Frequency (MHz)', fontsize=16)
plt.ylabel('Intensity (Kelvins)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
