#Python Script to deterine whether an event is long or short
#plot the raw data from the observation
#HISTORY
#19JUL15 GIL initial version based on e.py
#
import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
import numpy as np

dy = -1.

nargs = len( sys.argv)

if nargs < 4:
    print "G: Gauge event length.   If a time series with RMS > min and duration > N"
    print " Then this is a long event.   Count short and long events.  "
    print "Usage: "
    print "G N_Sigma  N_short <event files> "
    print " Where:"
    print " N_Sigma  Number of Sigma of samples to declare an Event"
    print " N_Short  Maximum number of Samples to consider a short event"
    exit()

nSigma = np.float(sys.argv[1])
nShort = np.int(sys.argv[2])
iFile = 3

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

#for symbol, value in locals().items():
#    print symbol, value

note = "Events"
nplot = 0
allNoiseCount = 0
longEventCount = 0
shortEventCount = 0

lastRmsEvent = 0
shortCountsPerHour = np.zeros(24)
longCountsPerHour = np.zeros(24)
allNoiseCountsPerHour = np.zeros(24)

# max number of high RMS intervals consider an event
maxShort = 10

for iii in range(iFile, min(nargs,40000)):

    filename = sys.argv[iii]

    rs = radioastronomy.Spectrum()
#    print filename
    rs.read_spec_ast( filename)
    rs.azel2radec()    # compute ra,dec from az,el

#    print("GAL Lon,Lat: %8.3f, %8.3f"  % (rs.gallon, rs.gallat))

    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    parts = aname.split('T')
    date  = parts[0]
    time  = parts[1]
    hourparts =  time.split('_')
    hourstr = hourparts[0]
    hourdigits = hourstr[0:2]
#    print "Time: Hours", time, hourparts, hourstr, hourdigits
    hour = int( hourdigits)
#    time  = time.replace('_',':')
    
    gallon = rs.gallon
    gallat = rs.gallat
    label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( time,rs.telaz,rs.telel,gallon,gallat)
    xs = rs.xdata 
    ya = rs.ydataA
    yb = rs.ydataB
    if rs.nTime < 1:
        print "Not an Event: ",filename
        continue

    xv = np.zeros(rs.nSamples*2)
    yv = np.zeros(rs.nSamples*2)
    j = 0
    dt = 0.5/rs.bandwidthHz
    t = -2. * dt * rs.refSample
    if iii == -1:  # no op
        print "First Time %12.9f (s); Delta T = %12.9f (s)" % (t, dt)

    for i in range(rs.nSamples):
        yv[j] = ya[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        yv[j] = yb[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)

    ymin = min(yv)
    ymax = max(yv)
    ymed = (ymin+ymax)*.5
    count = rs.count

# count number of possible short intervals in the time sequence
    nIntervals = int(2*rs.nSamples/nShort)
    shortCount = 0
    noEventCount = 0
    rmss = np.zeros( nIntervals)
    mins = np.zeros( nIntervals)
    maxs = np.zeros( nIntervals)
    j0 = 0
    jend = nShort
    rmsNoEvent = 0.
    rmsEvent = 0.

    for j in range( nIntervals):
#        aStd = np.sqrt(np.std( yv[j0:jend]))
        aStd = np.std( yv[j0:jend])
        mins[j] = -np.min( yv[j0:jend])
        maxs[j] = np.max( yv[j0:jend])
        rmss[j] = aStd
        j0 = j0 + nShort
        jend = jend + nShort

    # assume the median of events is a non event
    rmsMedian = np.median(rmss)
    eventSignal = nSigma*rmsMedian

    for j in range( nIntervals):
        aStd = rmss[j]
        if ((maxs[j] > eventSignal) or (mins[j] > eventSignal)):
            shortCount = shortCount + 1
            rmsEvent = rmsEvent + aStd
        else:
            noEventCount = noEventCount + 1
            rmsNoEvent = rmsNoEvent + aStd
            
    # compute index to middle event candidate, plus a couple samples
    jmiddle = int(nIntervals/2)+2
#    print "nIntervals: %d; Event Signal: %7.5f" % (nIntervals, eventSignal)
#    print "Mins: %7.5f %7.5f %7.5f" % (mins[2], mins[jmiddle], mins[nIntervals-2])
#    print "Maxs: %7.5f %7.5f %7.5f" % (maxs[2], maxs[jmiddle], maxs[nIntervals-2])

    if shortCount > 0:
        rmsEvent = rmsEvent/shortCount
    if noEventCount > 0:
        rmsNoEvent = rmsNoEvent/noEventCount

# if the same data was exactly the same as previousm, then the file was duplicated.
    if lastRmsEvent == rmsEvent:
        continue
    lastRmsEvent = rmsEvent   # remember last RMS

    t = shortCount * dt * nShort

    if rmsEvent < 4.*rmsNoEvent:
        allNoiseEvent = True
    else:
        allNoiseEvent = False

    if shortCount > maxShort or allNoiseEvent:
        shortEvent = False
#        print "%s Long  Event, Duration %9.6f (%4d, %7.3f, %7.3f)" % (filename, t, shortCount, rmsEvent, rmsNoEvent)
    else:
        print "%s Short Event, Duration %9.6f (%4d, %7.3f, %7.3f)" % (filename, t, shortCount, rmsEvent, rmsNoEvent)
        shortEvent = True
    
#    print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, rmsMedian, shortCount, label))
    if nplot <= 0:
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(date)
        nplot = nplot + 1
    note = rs.noteA
#    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
    plt.ylim(yallmin,1.25*yallmax)

    if shortEvent:
        shortEventCount = shortEventCount + 1
        if nplot < 10:
            plt.plot(xv, yv, colors[nplot-1], linestyle=linestyles[nplot-1],label=label)
            nplot = nplot + 1
        shortCountsPerHour[hour] = shortCountsPerHour[hour] + 1
    else:
        if allNoiseEvent:
            allNoiseCount = allNoiseCount + 1
            allNoiseCountsPerHour[hour] = allNoiseCountsPerHour[hour] + 1
        else:
            longEventCount = longEventCount + 1
            longCountsPerHour[hour] = longCountsPerHour[hour] + 1

plt.title(note)
plt.xlabel('Time Offset From Event (s)')
plt.ylabel('Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()

print "Short Event Count: %d; Long Event Count: %d; All Noise Count: %d" % ( shortEventCount, longEventCount, allNoiseCount)

print "Maximum Duration to consider a Short Event: %12.9f (s)" % (maxShort * dt * nShort)

# now report counts per hour
print "Summary of Types per hour:"
print ""
print "Hour Short  Long   All Noise"
for i in range(24):
    print " %2d  %4d  %4d    %4d" % ( i, shortCountsPerHour[i], longCountsPerHour[i], allNoiseCountsPerHour[i])
