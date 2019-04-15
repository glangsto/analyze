#Python Script to plot and the fourier transform of blocks of raw NSF events 
#HISTORY
#19APR11 GIL initial version for folding many files
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
ifile = 1
nfiles = nargs-ifile

MAXEVENTS = 5000

if nfiles < 2:
    print "MATCH: Match events listed in Two directories"
    print "Usage: MATCH dir1 dir2"
    exit()
    

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

    print "File  Index: ", ifile

    dir1 = sys.argv[ifile]
    dir2 = sys.argv[ifile+1]
    print "Dir 1: ", dir1
    print "Dir 2: ", dir2
    gridtype = 'PULSAR'

    from os import listdir
    from os.path import isfile, join
    files1 = [f for f in listdir(dir1) if isfile(join(dir1, f))]
    files2 = [f for f in listdir(dir2) if isfile(join(dir2, f))]

    count = 0
    print "%5d Files in Directory 1" % (len(files1))
#    print files
    print "%5d Files in Directory 2" % (len(files2))
    print files2
    nEve1 = 0
    nEve2 = 0

    rs = radioastronomy.Spectrum()

    event1s = files1
    event2s = files2

    for filename in files1:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "eve":
            print filename
            event1s[nEve1] = filenames
            nEve1 = nEve1 + 1
        else:
            continue

    utc1s = np.zeros(nEve1)

    iii = 0
    for filename in files1:
        fullname = join(dir1, filename)
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        utc1s[iii] = rs.utc
        iii = iii + 1


    for filename in files2:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "eve":
            print filename
            nEve2 = nEve2 + 1
        else:
            continue


        

if __name__ == "__main__":
    main()

