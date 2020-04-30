#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#20Apr01 GIL update for new python environment
#20FEB15 GIL normalize for different integration times
#19DEC30 GIL add title option
#17NOV21 GIL use time in file to show date and time
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import matplotlib.pyplot as plt
import sys
import numpy as np
import radioastronomy
import interpolate

dy = -1.

nargs = len( sys.argv)
verbose = True

# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False # flag interpolate over spike in center of spectrum
doDebug = False   # flag printing debug info
doSave = False    # flag saving intermediate files
flagRfi = True    # flag flagging RFI
doFold = False    # fold spectra to address an old issue; not normally used.
# put your list of known RFI features here.  Must have at least two, if flagRfi is true.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]   # integer number of channels to interpolate over
outDir = "./"     # define output directory for saving files
myTitle = ""      # Default no input title

# currently used velocities for baseline fitting
maxvel = 180.
minvel = -550.

nargs = len(sys.argv)
#first argument is the averaging time in seconds
namearg = 1
iarg = 1          # start searching for input flags
# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-F':
        print('Folding specectra')
        doFold = True
    elif sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-D':
        print('Adding Debug Printing')
        doDebug = True
    elif sys.argv[iarg].upper() == '-S':
        print('Saving Average Spectra in -ave files')
        doSave = True
    elif sys.argv[iarg].upper() == '-B':
        print('Baseline subtraction')
        doSub = True
    elif sys.argv[iarg].upper() == '-O':   # now look for flags with arguments
        iarg = iarg+1
        namearg = namearg+1
        outDir = sys.argv[iarg]
        print('Output directory: ', outDir)
    elif sys.argv[iarg].upper() == '-T':   # if plot title provided
        iarg = iarg+1
        namearg = namearg+1
        myTitle = sys.argv[iarg]
        print('Plot Title : ', myTitle)
    elif sys.argv[iarg].upper() == '-VA':   # now look for flags with arguments
        iarg = iarg+1
        namearg = namearg+1
        minvel = float(sys.argv[iarg])
        print('Minimum velocity for baseline fit: %7.2f km/sec ' % (minvel))
    elif sys.argv[iarg].upper() == '-VB':   # now look for flags with arguments
        iarg = iarg+1
        namearg = namearg+1
        maxvel = float(sys.argv[iarg])
        print('Maximum velocity for baseline fit: %7.2f km/sec ' % (maxvel))
    else:
        break
    iarg = iarg + 1
    namearg = namearg+1
# end of while not reading file names

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
scalefactor = 1.0
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9

linelist = [1420.0, 1418.0]  # RFI lines in MHz
linewidth = [7, 7]

#for symbol, value in locals().items():

# initialize spectrum for reading and plotting
rs = radioastronomy.Spectrum()

nplot = 0

# plot no more than N spectra
for iii in range(namearg, min(nargs,20)):

    filename = sys.argv[iii]
    if verbose:
        print('%5d: %s' % (iii, filename))

#    print(filename)
    rs.read_spec_ast( filename)
# for averages can not use az,el to get ra,dec and glat, glon
#    rs.azel2radec()    # compute ra,dec from az,el 

#    print("GAL Lon,Lat: %8.3f, %8.3f"  % (rs.gallon, rs.gallat))


    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
# now compute strings for plotting
    strtime = rs.utc.isoformat()
    parts = strtime.split('T')
    date  = parts[0]
    time  = parts[1]
    time  = time.replace('_',':')
    parts  = time.split('.')
    time = parts[0]
    
    gallon = rs.gallon
    gallat = rs.gallat
    label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( time,rs.telaz,rs.telel,gallon,gallat)
    xv = rs.xdata  * 1.E-6 # convert to MHz
    nData = len( xv)
    n6 = int(nData/6)
    n56 = 5*n6

    # normize for different integration times
    rs.ydataA = rs.ydataA/rs.count
    yv = rs.ydataA

    # The latest versions of the software had a different normalization
    ymedian = np.median(yv[n6:n56])
    # if the latest software, the scale factor is just 1.
    if ymedian > .001:
        scalefactor = 1.0

    yv = rs.ydataA * scalefactor
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)

    yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

    ymin = min(yv)
    ymax = max(yv)
    ymed = np.median(yv)
    count = rs.count

    print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
    if nplot <= 0:
        fig,ax1 = plt.subplots(figsize=(10,6))
#        plt.hold(True)
        fig.canvas.set_window_title(date)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 

        nplot = nplot + 1
#    note = rs.noteA
    note = rs.site
#    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
# scale min and max intensities for nice plotting
    plt.ylim(0.9*yallmin,1.25*yallmax)

    plt.plot(xv, yv, colors[iii-1], linestyle=linestyles[iii-1],label=label, lw=2)
if myTitle == "":
    myTitle = note
plt.title(myTitle, fontsize=16)
plt.xlabel('Frequency (MHz)',fontsize=16)
ylabel = 'Intensity (%s)' % rs.bunit
plt.ylabel(ylabel, fontsize=16)
plt.legend(loc='upper right')
plt.show()
