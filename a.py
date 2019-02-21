#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class

#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import matplotlib.pyplot as plt
import sys
import datetime
import statistics
import radioastronomy

avetimesec = 1800.
dy = -1.

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

#for symbol, value in locals().items():
#    print symbol, value

nplot = 0
nhot = 0
ncold = 0

# first read through all data and find hot load
for iii in range(1, nargs):

    filename = sys.argv[iii]

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast( filename)
    rs.azel2radec()    # compute ra,dec from az,el

    if rs.telel < 0:
        if nhot == 0:
            hot = rs
            nhot = 1
        else:
            hot.ydataA = hot.ydataA + rs.ydataA
            hot.count = hot.count + rs.count
            nhot = nhot + 1

if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA / float(nhot)
    print "Found %3d hot load obs" % nhot
else:
    print "No hot load data, can not calibrate"
    exit()

xv = hot.xdata * 1.e-6 # convert to MHz
yv = hot.ydataA
yallmin = min(yv)
yallmax = max(yv)

fig,ax1 = plt.subplots(figsize=(10,6))
plt.hold(True)
fig.canvas.set_window_title(filename)
az = hot.telaz
el = hot.telel
label = 'Hot Load Average'
ymin = min(yv)
ymax = max(yv)
ymed = statistics.median(yv)
count = hot.count
ncold = 0
print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
plt.plot(xv, yv, colors[0], linestyle=linestyles[0],label=label)

avetime = datetime.timedelta(seconds=avetimesec)

nread = 0        
# now read through all data and average cold sky obs
for iii in range(1, nargs):

    filename = str(sys.argv[iii])

    rs = radioastronomy.Spectrum()
#    print filename
    rs.read_spec_ast( filename)
    rs.azel2radec()    # compute ra,dec from az,el

    if rs.telel > 0:
        if ncold == 0:
            cold = rs
            ncold = 1
            # sums are weighted by durations
            cold.ydataA = cold.ydataA * cold.durationSec
            # keep track of observing time for weighted sum
            timesum = cold.durationSec
        else:
            # time difference is between mid-points of integrations
            dt = rs.utc - cold.utc 
            # add the time since midpoint of latests
            dt = dt + datetime.timedelta(seconds=rs.durationSec/2.)
            # plus time before start of the first
            dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)
            # if time to average
            if dt > avetime:
                cold.ydataA = cold.ydataA/float(timesum)

                gallon = cold.gallon
                gallat = cold.gallat
                az = cold.telaz
                el = cold.telel
                parts = filename.split('/')
                nparts = len(parts)
                aname = parts[nparts-1]
                parts = aname.split('.')
                aname = parts[0]
                parts = filename.split('/')
                nparts = len(parts)
                aname = parts[nparts-1]
                parts = aname.split('.')
                aname = parts[0]
                parts = aname.split('T')
                date  = parts[0]
                time  = parts[1]
                time  = time.replace('_',':')

                label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % (time, az,el,gallon,gallat)
                xv = cold.xdata * 1.e-6
                yv = cold.ydataA * scalefactor

                xmin = min(xv)
                xmax = max(xv)
                xallmin = min(xmin,xallmin)
                xallmax = max(xmax,xallmax)
                ymin = min(yv)
                ymax = max(yv)
                ymed = statistics.median(yv)
                count = cold.count
                ncold = 0
                print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
                nplot = nplot + 1
                note = cold.noteA
                    #print('%s' % note)
                yallmin = min(ymin,yallmin)
                yallmax = max(ymax,yallmax)
                ncolor = min(nmax-1, nplot) 
                plt.plot(xv, yv, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            else: # else ont enough time yet, average cold data
                cold.count = cold.count + rs.count 
                cold.ydataA = cold.ydataA + (rs.ydataA * cold.durationSec)
            # keep track of observing time for weighted sum
                timesum = timesum + cold.durationSec
            # end if not a enough time
        # end if a cold file
    #end for all files to sum

plt.xlim(xallmin,xallmax)
plt.ylim(0,1.5*yallmax)
plt.title(date)
fig.canvas.set_window_title(date)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
plt.xlabel('Frequency (MHz)', fontsize=14)
plt.ylabel('Intensity (Counts)', fontsize=14)
plt.legend(loc='upper right')
plt.show()
