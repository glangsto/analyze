#Python Script to plot the fourier transform of blocks of raw NSF events 
#plot the FFT of raw data from an event
#HISTORY
#19FEB21 GIL initial version
#
import matplotlib.pyplot as plt
import sys
import statistics
import radioastronomy
import interpolate

dy = -1.

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

# separate arguments form file names
ifile = 1
iii = ifile

nblock = 128
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-BL":
        nblock = np.int( sys.argv[iii+1])
        iii = iii + 1
        print "FFT Block Size: ", nblock
        aFix = True
        ifile = ifile + 2
    iii = iii + 1

note = "Events"
nplot = 0
nfiles = nargs-ifile
for iii in range(1, min(nfiles,25)):

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
    time  = time.replace('_',':')
    
    gallon = rs.gallon
    gallat = rs.gallat
    label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( time,rs.telaz,rs.telel,gallon,gallat)
    xs = rs.xdata 
    if rs.nTime < 1:
        print "Not an Event: ",filename
        continue

    xv = np.zeros(rs.nSamples*2)
    yv = np.zeros(rs.nSamples*2)
    nSamples = 2L * rs.nSamples
    j = 0
    dt = 0.5/rs.bandwidthHz
    t = xs[0]
    for i in range(rs.nSamples):
        yv[j] = rs.ydataA[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        yv[j] = rs.ydataB[i]
        xv[j] = t
        j = j + 1
        t = t + dt

    # will compute several FFTs, discarding any extra samples at end
    nffts = nSamples/nblock
    lll = 0
    for jjj in range(nffts):
        for kkk in range(nblock):
            ablock[kkk] = yv[lll]
            lll += 1

    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)

    ymin = min(yv)
    ymax = max(yv)
    ymed = statistics.median(yv)
    count = rs.count

    print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
    if nplot <= 0:
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(date)
        nplot = nplot + 1
    note = rs.noteA
#    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
#    plt.ylim(0.9*ymin,1.5*yallmax)
#    plt.ylim(0.9*ymin,1.25*yallmax)
    plt.ylim(0.,1.25*yallmax)


    plt.plot(xv, yv, colors[iii-1], linestyle=linestyles[iii-1],label=label)
plt.title(note)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
