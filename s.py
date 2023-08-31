#Python Script to plot summarize NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#23AUG16 GIL check for big changes to Tsys 
#19SEP11 GIL no longer use statistics
#18MAR05 GIL only summarize .hot and .ast files
#16SEP30 GIL show frequency and bandwidth
#16AUG29 GIL make more efficient
#
import sys
import numpy as np
import radioastronomy
import gainfactor as gf

percent = 10           # default percentage change to flag
nargs = len( sys.argv)
if nargs < 2:
    print('S: Summarize a list of observations')
    print('S usage:')
    print('S [-s %] [-v] <filenames>')
    print('where')
    print('-v          optionally print verbose summary information')
    print('-s %        optionally report jumps than % (range 1 to 100)')
    print('            Default percentage: %7.1f %% ' % (percent))   
    print('<filenames> list of file names to summarize')
    print('  Only *.ast, *.hot and *.cld files will be read. Others will be skipped')
    exit()

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

#for symbol, value in locals().items():
#    print symbol, value
lastaz = 0.
lastel = 0.
lastfreq = 0.
lastbw = 0
lastgain = 0.
nsame = 0
nRead = 0

labelfmt = '%s %5s,%5s %5.1f,%5.1f: %8.2f,%5.2f %5.1f - %s' 
names = sys.argv[1:nargs]
names = sorted(names)

##label = labelfmt % (time,rs.telaz,rs.telel,gallon,gallat,freq,bw,gain, filename)
rs = radioastronomy.Spectrum()

# check inputs twice to allow optional inputs in either order
if names[0] == '-s':
    percent = float(names[1])
    names = names[2:]
                    
if names[0] == '-v':
    verbose = True
    names = names[1:]
else:
    verbose = False

if names[0] == '-s':
    percent = float(names[1])
    names = names[2:]
                    
if names[0] == '-v':
    verbose = True
    names = names[1:]

# prepare to store Tsys
nNames = len(names)
names2 = names    # 2nd array to hold in case any files are skipped
tSys = np.zeros(nNames)
tRms = np.zeros(nNames)
velMin = 200.*1000.   # m/sec
velMax = 300.*1000   #  m/sec

iMin = 0      # index for computing tSys, tRms
iMax = 0      #

# for all input names
for filename in names:

    # split out directory names
    parts = filename.split('/')
    nparts = len(parts)
    # last part of name is the file name
    aname = parts[nparts-1]
    parts = aname.split('.')
    # determine extension type
    nparts = len(parts)
    # if no extension, 
    if nparts < 2: 
        continue  ### then not a file to summarize
    
# only summarize astronomy files
    if (parts[1] != 'ast') and (parts[1] != 'hot') and (parts[1] != 'cld') \
       and (parts[1] != 'jy') and (parts[1] != 'kel'):
        continue
    if verbose:
        print('File: ',filename)
    rs.read_spec_ast( filename)
    rs.azel2radec()    # compute ra,dec from az,el
    gallon = rs.gallon
    gallat = rs.gallat
    freq = rs.centerFreqHz / 1.E6
    bw = rs.bandwidthHz / 1.E6
    gain = rs.gains[0]
    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    adatetime = str(rs.utc)
    parts = adatetime.split(' ')
    date  = parts[0]
    time  = parts[1]
    parts  = time.split('.')
    time = parts[0]

    # now compute the indices for the velocity range for RMS and TSYS
    freqHz = np.array(rs.xdata)
    if verbose:
        print("freqHz: %.1f  (%d)" % (freqHz[0], len(freqHz)))
    velmps = gf.velocity( freqHz, rs.refFreqHz)
    
    # next comput indicies for tsys, trms cal
    names2[nRead] = filename
    ixMin, ixMax = gf.velocity_to_indicies( velmps, velMin, velMax)
    tSys[nRead] = np.median(rs.ydataA[ixMin:ixMax])
    tRms[nRead] = np.std(rs.ydataA[ixMin:ixMax])
    
    if rs.telaz != lastaz or rs.telel != lastel or lastbw != bw or lastfreq != freq or lastgain != gain:
        lastbw = bw
        lastfreq = freq
        lastaz = rs.telaz
        lastel = rs.telel
        lastgain = gain

        if nRead == 0:
            label = (labelfmt % (time,rs.telaz,rs.telel,gallon,gallat,freq,bw,gain,filename))
        elif nsame > 1:
            print("%4d %s " % (nsame, label))
        label = labelfmt % (time,rs.telaz,rs.telel,gallon,gallat,freq,bw,gain, filename)
        if nRead == 0:
            print('Count  Time    Az    El   G-Lon G-Lat  Frequency  BW   Gain    Filename')
        print("%4d %s " % (1, label))
        nsame = 1
    else:
        nsame = nsame + 1
    nRead = nRead + 1

    label = labelfmt % (time,rs.telaz,rs.telel,gallon,gallat,freq,bw,gain, filename)

if nsame > 0:
    print("%4d %s " % (nsame, label))

print("")
# now find jumps in intensity
iWidth = 3
if (3*iWidth) > nRead:   # if not enough data to look for jumps
    exit()

fraction = percent/100.

print('Looking for jumps in intensity > %7.1f %%' % (percent))

for iii in range(iWidth,(nRead-iWidth)):
    # get the temperatures in counts
    t1 = np.median(tSys[(iii-iWidth):iii])
    dt1 = np.std(tSys[(iii-iWidth):iii])
    t2 = np.median(tSys[iii:(iii+iWidth)])
    dt2 = np.std(tSys[iii:(iii+iWidth)])

    if dt1 < 1.:
        dt1 = 1.   # minimum noise is 1 count
    if dt2 < 1.:
        dt2 = 1.   # minimum noise is 1 count


    # offset counts between the two sides
    dt = np.abs(t1-t2)
    if (dt > (fraction*t1)) or (dt > (fraction*t2)):
        print("At %s, Count Jump; %7.2f+/-%5.2f =! %7.2f+/-%5.2f" % \
              (names2[iii], t1, dt1, t2, dt2))
    
