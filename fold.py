#Python Script to plot and the fourier transform of blocks of raw NSF events 
#HISTORY
#19APR12 GIL fix average Time calculation
#19APR11 GIL improve labeling
#19APR10 GIL initial version based on FFT
#
import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
from scipy.fftpack import fft
import numpy as np
import copy
import GridClass
from scipy.signal import blackman

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
ncolor = min( len(colors), 7)

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

# separate arguments form file names
iCount = 1
ifile = 2
iii = ifile

nday = 24             # by default divide day in 24 hours
nblock = 256          # number of samples to FFT
sigma = 5.0           #
kpercount = 1.0       # calibration into Kelvin units
note = ""             # optional note for top of plot

# read through arguments extracting parameters
while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-BL":
        nblock = np.int( sys.argv[iii+1])
        iii = iii + 1
        print "FFT Block Size: ", nblock
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-KP":
        kpercount = np.float( sys.argv[iii+1])
        iii = iii + 1
        print "Kelvins Per Count: ", kpercount
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-ND":
        ndays = np.int( sys.argv[iii+1])
        iii = iii + 1
        print "Divide Day into N Parts: ", nday
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:3]) == "-SI":
        sigma = np.float( sys.argv[iii+1])
        iii = iii + 1
        print "Keeping Events > %7.2f Sigma " % (sigma)
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print "Note: ", note
        ifile = ifile + 2
        aFix = True
    iii = iii + 1

iCount = ifile - 1

N = nblock                 # abreviation
ablock = np.zeros(nblock)  # creat array to FFT
eventblock = np.zeros(nblock)  # creat array to FFT
eventCounts = np.zeros(nday)  # count events per fraction of a day
eventAveGLat = np.zeros(nday) # count events per fraction of a day
eventAveGLon = np.zeros(nday) # count events per half hour
eventAveRa = np.zeros(nday) # count events per fraction of a day
eventAveDec = np.zeros(nday) # count events per half hour
w = blackman(N)
nu = np.zeros(nblock)  # creat frequency array
nplot = 0
ifile = 2
nfiles = nargs-ifile

print "Number of Files:                ",nfiles
if nfiles < 1:
    print "FOLD: Fourier Transform and sum a time series of events"
    print "Usage: FOLD <n channels] [-note <note for plot title>] <file 1>"
    print "Where optionally the following paramters may be applied"
    print " <n samples>     - Number of samples to FFT to produce a spectra (power of 2)"
    print " -kpercount <factor> - Gain factor to convert counts to Kelvin"
    print "   Estimate by an assumed system temperature for the band pass"
    print "   And running FFT without this factor to get the value in counts"
    print " -note <text>        - Note for the top of the plot"
    exit()
    
files = sys.argv[ifile:]
print "First File     :                ",files[0]


def gridAFile( nChan = 128, filename = ""):
    """
    Compute a grided image of the fourier transform of a time series
    Inputs: filename - name of the Event file to read
    nblock - number of samples to use to compute a single Fourier Transform
    """

    maxMagnitude = 0.
    maxEvent = radioastronomy.Spectrum()
    maxFile = ""

    nblock = 2 * int(nChan)
    nblock1 = nblock - 1
    N = nblock                 # abreviation
    N2 = N//2
#    print "N, N2: ", N, N2

    yp = np.zeros(N2)
    yp2 = np.zeros(N2)
    ysum = np.zeros(N2 - 1)
    nsum = 0      # count total number of spectra summed

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
#    time  = time.replace('_',':')
    
    gallon = rs.gallon
    gallat = rs.gallat
    label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( time,rs.telaz,rs.telel,gallon,gallat)
    xs = rs.xdata 
    if rs.nTime < 1:
        print "Not an Event: ",filename
        return 0

    # now reorganize data into a single time series
    xv = np.zeros(rs.nSamples*2)
    yv = np.zeros(rs.nSamples*2)
    yc = np.zeros(rs.nSamples,dtype=np.complex)
    ymag = np.zeros(rs.nSamples)
    nSamples = 2L * rs.nSamples
    j = 0
    dt = 0.5/rs.bandwidthHz
    t = xs[0]
    for i in range(rs.nSamples):
        yv[j] = rs.ydataA[i]
        yc[i] = rs.ydataA[i] + 1j*rs.ydataB[i]
        xv[j] = t
        j = j + 1
        t = t + dt
        yv[j] = rs.ydataB[i]
        xv[j] = t
        j = j + 1
        t = t + dt

    # compute magnitude of I/Q samples with vector math
    ymag = np.absolute(yc)
    # now find the maximum event in data series
    ymagmax = max(ymag)
    # compute RMS of entire time series
    yrms = np.sqrt(yv.dot(yv)/yv.size)
    # if a valid, noisy observation
    if yrms > 0. :
        nsigma = ymagmax/yrms
        if nsigma > maxMagnitude:
            maxEvent = copy.deepcopy(rs)
            maxMagnitude = nsigma
            maxFile = filename
    else:
        print "Problem observation: %7.3f rms in file %s\n", yrms, filename
    # end else not a signficant event

    # will compute several FFTs, discarding any extra samples at end
    nfft = int(nSamples/nblock) - 1
    lll = 0
# instead start at half a block to center event
    lll = int(nblock/2)
    # 
    BW = 1.E-6*rs.bandwidthHz
    # average time for spectra without event
    aveTime = (nfft-1)*nblock*1.e-6/(2.*BW)
    # frequency axis is always the same
    nu = np.linspace(0.0, BW, nblock/2) + (1.E-6*(rs.centerFreqHz-rs.bandwidthHz/2.))
    xp = nu[1:N2]
    xmin = min(nu)
    xmax = max(nu)
    # zero sum of all spectra in this file


    #create the grid with map parameters
    mygrid = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=0, ymax=nfft, width=nChan, \
                                height=nfft, dpi=1, FWHM=1., \
                                projection="Mercator", gridtype="")

    for jjj in range(nfft):
        J = float(jjj)
        # for each block of samples
        for kkk in range(nblock):
            ablock[kkk] = yv[lll]
            lll += 1
        yf = fft(ablock*w)  # apply blackman window/weighting to samples
#        yf = fft(ablock)   # fft without window
        yp = 2.0/N*np.abs(yf[1:N2])
        
        ymin = min(yp)
        ymax = max(yp)

        # unwrap spectra
        nnn = ((N2/2) - 1)
        for kkk in range(N2-1):
#            if jjj == nfft/2:
#                print "kkk, nnn: ",kkk,nnn
            yp2[kkk] = yp[nnn]
            nnn -= 1
            if nnn < 0:
                nnn = N2 - 2

        for mmm in range( N2-1):
            mygrid.convolve( xp[mmm], J, yp2[mmm], 1.0)

        if jjj == 0:
            print(' Max: %9.3f Min: %9.3f ; %3d %s' % (ymax*kpercount, ymax*kpercount, nfft, label))
            nyp = len(yp2)
            ysum = np.zeros(nyp)

# plot the middle spectrum for each events, until out of colors
        if jjj == (nfft/2):
            nplot = 1
            yp2 = kpercount * yp2
            yp3 = yp2[1:N2]
            label = '%s %s Lon,Lat=%5.1f,%5.1f' % ( date, time, gallon, gallat)
            plt.plot(xp, yp3, colors[nplot-1], linestyle=linestyles[nplot-1],label=label)
            note = rs.noteA
            print "Number of FFTs per time series: %5d" % ( nfft)
            print "Duration of Observations      : %12.6f s" % ( aveTime)
        else:
#            print "Ysum, yp2: ", len(ysum), len(yp2)
            ysum = ysum + yp2
            nsum += 1

    # now average of all spectra not during the peak event
    ysum = kpercount * ysum / float(nsum)
    ysum = ysum[1:N2] 
    avelabel = "Average (%8.6f s)" % aveTime
    plt.plot(xp, ysum, colors[1], linestyle=linestyles[1],label=avelabel)

# put the average in the top row of the image
    J = float(nfft)
    for mmm in range( N2-1):
        mygrid.convolve( xp[mmm], J, ysum[mmm], 1.0)
    print "Event Time         : ", maxEvent.utc
    print "Event File         : ", maxFile
    print "Maximum Event Sigma: %8.2f" % ( maxMagnitude)
    print "Event RA           :   %8.4f d Dec %8.4f d" % (maxEvent.ra, maxEvent.dec)
    print "Event G Lon        :   %6.2f d  Lat %6.2f d" % (maxEvent.gallon, maxEvent.gallat)
    
    datetime = "%s" % ( maxEvent.utc)
    parts = datetime.split('.')
    date = parts[0]

#        plt.xlim(xallmin,xallmax)
    plt.title(note)
    plt.xlabel('Frequency (MHz)')
    if kpercount == 1.:
        plt.ylabel('Intensity (Counts)')
        #            plt.ylim(yallmin,1.25*yallmax)
    else:
        plt.ylabel('Intensity (Kelvins)')
        #            plt.ylim(kpercount*yallmin,kpercount*1.25*yallmax)
    plt.legend(loc='upper right')
    plt.show()

    return nChan, nfft, xmin, xmax, mygrid

def main():
    """
    Main executable for gridding astronomical data
    """
    dpi = 1

    nargs = len(sys.argv)
    if nargs < 2:
        print 'fold : Fold A time series'
        print 'usage: fold <nChannels>  <eventfile>'
        exit()

    print "Count Index: ", iCount
    print "File  Index: ", ifile

    nChan = int(sys.argv[iCount])
    names = sys.argv[ifile:]
    print "First Name: ", names[0]
    gridtype = 'PULSAR'

    #create the grid with map parameters
    mywidth, myheight, xmin, xmax, mygrid = gridAFile( nChan, names[0])

    count = 0

    mygrid.normalize()
#    mygrid.check()
    zmin = 00.
    zmax = 2000.
# limit grid intensities for plotting
#    mygrid.limit(zmin, zmax)

    ticks = np.arange(0, mywidth, 30*dpi)
    ymin = 0
    ymax = myheight

    # round to compute number of MHz on plot
    nMHz = int((xmax - xmin)+.01)
    if nMHz < 3:   # If not many lablels 
        nMHz = 4
    elif nMHz > 8: # else too many labels
        nMHz = 8
    
    # these ticks are in "pixel" units
    xticks = np.arange(0, mywidth, int(mywidth/nMHz))

    # these ticks are in astronomical units (ie MHz)
    x_ticks = xmin + ((xticks*(xmax-xmin))/mywidth)
    nTicks = len(x_ticks)
    # ticks are on the even 10ths of MHz
    for iii in range(nTicks):
        x_ticks[iii] = int(10*(x_ticks[iii]+.095))/10.

    plt.imshow(mygrid.image, interpolation='nearest', cmap=plt.get_cmap('jet'))

    plt.title("Pulsar-Folded Time Series of Spectra")

    if gridtype == 'PULSAR':
        plt.xlabel("Frequency (MHz)")
        plt.ylabel("Spectrum Number")
        yticks = np.arange(0, myheight, 15*dpi)
        y_ticks = ymax - (ymax-ymin)*yticks/myheight
# diagnostics to figure out labeling
#        print "Y"
#        print yticks, y_ticks
#        print "X"
#        print xticks, x_ticks

        plt.yticks(yticks, y_ticks)
        plt.xticks(xticks, x_ticks, rotation='horizontal')
        plt.colorbar()

    plt.show()

if __name__ == "__main__":
    main()

