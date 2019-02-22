#Python Script to plot the fourier transform of blocks of raw NSF events 
#plot the FFT of raw data from an event
#HISTORY
#19FEB21 GIL initial version
#
import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
from scipy.fftpack import fft
import numpy as np

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
note = ""
# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-BL":
        nblock = np.int( sys.argv[iii+1])
        iii = iii + 1
        print "FFT Block Size: ", nblock
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print "Note: ", note
        ifile = ifile + 2
        aFix = True
    iii = iii + 1

N = nblock                 # abreviation
ablock = np.zeros(nblock)  # creat array to FFT
nu = np.zeros(nblock)  # creat frequency array
ysum = np.zeros((nblock/2)-1)
nplot = 0
nfiles = nargs-ifile+1

for iii in range(1, min(nfiles,25)):

    filename = sys.argv[iii]

    rs = radioastronomy.Spectrum()
#    print filename
    rs.read_spec_ast( filename)
    if note != "":
        rs.noteA = note
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
    # 
    BW = 1.E-6*rs.bandwidthHz
    # frequency axis is always the same
    nu = np.linspace(0.0, BW, nblock/2) + (1.E-6*rs.centerFreqHz)
    xmin = min(nu)
    xmax = max(nu)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)

    for jjj in range(nffts):
        # for each block of samples
        for kkk in range(nblock):
            ablock[kkk] = yv[lll]
            lll += 1
        yf = fft(ablock)
        xp = nu[1:N//2]
        yp = 2.0/N*np.abs(yf[1:N//2])
        ysum = ysum + yp

        ymin = min(yp)
        ymax = max(yp)

        print(' Max: %9.2f Min: %9.2f ; %3d %s' % (ymax, ymax, jjj, label))
        plabel = ("%d" % (jjj))
        if nplot <= 0:
            fig,ax1 = plt.subplots(figsize=(10,6))
            fig.canvas.set_window_title(date)
        nplot = nplot + 1
#    print('%s' % note)
        yallmin = min(ymin,yallmin)
        yallmax = max(ymax,yallmax)

        plt.plot(xp, yp, colors[nplot-1], linestyle=linestyles[nplot-1],label=plabel)
note = rs.noteA
ysum = ysum / np.float(nffts)
plt.plot(xp, ysum, colors[0], linestyle=linestyles[0],label="Sum", lw=3)
plt.xlim(xallmin,xallmax)
plt.ylim(yallmin,1.25*yallmax)
plt.title(note)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Intensity (Counts)')
plt.legend(loc='upper right')
plt.show()
