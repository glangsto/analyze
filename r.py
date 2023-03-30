#Python Script to plot raw NSF record data.
#plot the raw data from the observation
#HISTORY
#23Mar04 GIL fix matplot lib changes/deprications.   Short label if el < 0
#22May20 GIL if file not recognized skip
#20Dec24 GIL remove 20 from plot file names
#20Aug28 GIL update plotting to a file
#20Jun16 GIL add plotting to a file
#20Apr01 GIL update for new python environment
#20FEB15 GIL normalize for different integration times
#19DEC30 GIL add title option
#17NOV21 GIL use time in file to show date and time
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import matplotlib as mpl
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
doPlotFile = False
plotFileDir = "~/"
# put your list of known RFI features here.  Must have at least two, if flagRfi is true.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
doFold = False
# fold spectra to address an old issue; not normally used.

# put your list of known RFI features here.  Must have at least two, if flagRfi is true.

plotFrequency = True
# define reference frequency for velocities (MHz)
nuh1 = 1420.40575 # neutral hydrogen frequency (MHz)
nuoh1= 1612.231   # OH line
nuoh2= 1665.402   # OH line
nuoh3= 1667.359   # OH Line
nuoh4= 1720.530   # OH Line

# select the frequency for plotting velocities
nuRefFreq = nuh1

linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]   # integer number of channels to interpolate over
outDir = "./"     # define output directory for saving files
myTitle = ""      # Default no input title
note = ""

# currently used velocities for plotting range
maxvel = 180.
minvel = -180.
maxPlot = int(25)
fileTag = ""
firstdate = ""
doScaleAve = False

nargs = len(sys.argv)
#first argument is the averaging time in seconds
namearg = 1
iarg = 1          # start searching for input flags
if nargs < 2:
    print("R: Plot Raw counts of telescope observations")
    print("Usage: R <flags> <files>")
    print("Where <flags> are:")
    print("-A optionally scale intensities by count of spectra averaged")
    print("-B <sample> Set first sample to plot (default is 1/4 of samples)")
    print("-C optionally flag the center of the band")
    print("-E <sample> Set last sample to plot (default is end of samples)")
    print("-H optionally set the high velocity region for baseline fit")
    print("-K optionall save average hot and cold load calibration observations")
    print("-L optionally set the low velocity region for baseline fit")
    print("-N <number> optionally set the number of spectra to plot")
    print("-P write PNG and PDF files instead of showing plot")
    print("-Q optionally plot intensity versus freQuency, instead of velocity")
    print("-S <filename> optionally set summary file name")
    print("-U optionally update reference frequency for a different line")
    print("   ie -U 1612.231, 1665.402, 1667.349, 1720.530 or 1420.40575")
    print("-V optionally plot velocity")
    print("-Z <file tag> optionally add tag to PDF and PNG file names")
    print("-F optionally fold spectra (obsolete)")
    print("")
    print("Glen Langston - NSF Dec 24, 2020")
    exit()

xa = -1
xb = -1

namearg = 1    
# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-F':
        print('Folding specectra')
        doFold = True
    elif sys.argv[iarg].upper() == '-A':
        toScalAve = True
    elif sys.argv[iarg].upper() == '-B':   # if setting beginning sample
        iarg = iarg + 1
        xa = int( sys.argv[iarg])
        print('Plotting starting at channel: %4d' % (xa))
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-E':   # if setting ending sample
        iarg = iarg + 1
        xb = int( sys.argv[iarg])
    elif sys.argv[iarg].upper() == '-H':
        iarg = iarg+1
        maxvel = np.float( sys.argv[iarg])
        print('Maximum (high) velocity for sum: %7.2f km/sec' % (maxvel))
    elif sys.argv[iarg].upper() == '-L':
        iarg = iarg+1
        minvel = np.float( sys.argv[iarg])
        print('Minium (low)  velocity for sum: %7.2f km/sec' % (minvel))
    elif sys.argv[iarg].upper() == '-N':   # if number of spectra to plot
        iarg = iarg+1
        maxPlot = int(sys.argv[iarg])
        if maxPlot < 1:
            print("Not Plotting")
        else:
            print("Plot will have a maximum of %d spectra" % (maxPlot))
    elif sys.argv[iarg].upper() == '-P':
        doPlotFile = True
        iarg = iarg+1
        plotFileDir = sys.argv[iarg]
    elif sys.argv[iarg].upper() == '-Q':
        plotFrequency = True
    elif sys.argv[iarg].upper() == '-R':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-T':   # if plot title provided
        iarg = iarg+1
        myTitle = sys.argv[iarg]
        print('Plot Title : ', myTitle)
    elif sys.argv[iarg].upper() == '-V':   # default is plotting Frequency
        plotFrequency = False              # plot velocity
    elif sys.argv[iarg].upper() == '-VA':   # now look for flags with arguments
        iarg = iarg+1
        minvel = float(sys.argv[iarg])
        print('Minimum velocity for baseline fit: %7.2f km/sec ' % (minvel))
    elif sys.argv[iarg].upper() == '-VB':   # now look for flags with arguments
        iarg = iarg+1
        maxvel = float(sys.argv[iarg])
        print('Maximum velocity for baseline fit: %7.2f km/sec ' % (maxvel))
    elif sys.argv[iarg].upper() == '-CEN':   # if nU ref is provided (in MHz)n
        iarg = iarg+1
        nuRefFreq = float(sys.argv[iarg])
        print( 'Reference Frequency : %9.3f MHz' % (nuRefFreq))
    elif sys.argv[iarg].upper() == '-Z':     # label written files
        iarg = iarg+1
        fileTag = str(sys.argv[iarg])
        print( 'File tag: %s' % (fileTag))
    else:
        break
#    namearg = iarg + 1
    namearg = iarg + 1
    iarg = iarg + 1
# end of while not reading file names

# to create plots in cronjobs, must use a different backend
if doPlotFile:
    mpl.use('Agg')
import matplotlib.pyplot as plt
import gainfactor as gf

if plotFrequency:
    print( "Ploting Intensity versus Frequency")
else:
    print( "Ploting Intensity versus Velocity")

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['g', 'b', 'r', 'b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g']
colors =  ['b','r','g', 'b','r','g','b','r','g','c','m','y','c','m','y','c','m','y','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g','b','r','g']
#colors = ['g', 'b', 'r', '-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

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

c = 299792.458  # (Speed of light  km/sec)
nplot = 0

# plot no more than N spectra
for iii in range(namearg, min(nargs,30)):

    filename = sys.argv[iii]
    if verbose:
        print('%5d: %s' % (iii, filename))

    rs.read_spec_ast( filename)
# for averages can not use az,el to get ra,dec and glat, glon
#    rs.azel2radec()    # compute ra,dec from az,el 

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

    if rs.nChan < 1:
        print("Skipping file: %s" % (filename))
        continue
    gallon = rs.gallon
    gallat = rs.gallat
    if rs.telel > 0:
        label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( time,rs.telaz,rs.telel,gallon,gallat)
    else:
        label = '%s, AZ,EL: %5s,%5s' % ( time,rs.telaz,rs.telel)
    xv = rs.xdata  * 1.E-6 # convert to MHz
    nData = len( xv)
    n6 = int(nData/6)
    n56 = 5*n6

    if xa < 0:
        xa = 0
    if xb < 0:
        xb = nData
        
    if firstdate == "":
        firstdate = date
    lastdate = date
    
    vel = np.zeros(nData)
    for jjj in range (0, nData):
        vel[jjj] = c * (nuRefFreq - xv[jjj])/nuRefFreq

    if not plotFrequency:
        xa, xb = gf.velocity_to_indicies( vel, minvel, maxvel)
    
    # normize for different integration times
    if doScaleAve:
        rs.ydataA = rs.ydataA/rs.count
    else:
        rs.ydataA = rs.ydataA/rs.nave
    yv = rs.ydataA

    # The latest versions of the software had a different normalization
    ymedian = np.median(yv[n6:n56])
    # if the latest software, the scale factor is just 1.
    if ymedian > .001:
        scalefactor = 1.0

    yv = rs.ydataA * scalefactor
    if not plotFrequency:
        xv = vel
    xmin = min(xv[xa:xb])
    xmax = max(xv[xa:xb])
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)

    if flagRfi:
        yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

    if flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
        icenter = int(nData/2)
        yv[icenter] = (yv[icenter-2] + yv[icenter+2])*.5
        yv[icenter-1] = (3.*yv[icenter-2] + yv[icenter+2])*.25
        yv[icenter+1] = (yv[icenter-2] + 3.*yv[icenter+2])*.25

    ymin = min(yv)
    ymax = max(yv)
    ymed = np.median(yv)
    count = rs.count

    print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label))
    if nplot <= 0 and maxPlot > 0:
#        fig = plt.figure(date)
        fig,ax1 = plt.subplots(figsize=(10,6))
#        fig.canvas.set_window_title(date)
#        plt.manager.canvas.set_window_title(date)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 

    nplot = nplot + 1
    if nplot > maxPlot:
        break
    note = rs.site
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
# scale min and max intensities for nice plotting
    plt.ylim(0.9*yallmin,1.25*yallmax)

    if plotFrequency:
        plt.plot(xv[xa:xb], yv[xa:xb], colors[nplot], linestyle=linestyles[iii-1],label=label, lw=2)
    else:
        plt.plot(xv[xa:xb], yv[xa:xb], colors[nplot], linestyle=linestyles[iii-1],label=label, lw=2)

# end for all names loop
if (maxPlot < 1) or (nplot < 1):
    print("No Plots, exiting")
    exit()

if myTitle == "":
    myTitle = note
    
plt.title(myTitle, fontsize=16)
plt.xlabel('Frequency (MHz)',fontsize=16)
ylabel = 'Intensity (%s)' % rs.bunit
plt.ylabel(ylabel, fontsize=16)
plt.legend(loc='upper right')
# if writing files
if doPlotFile:
    firstdate = firstdate[2:]
    if fileTag == "":
        fileTag = "R-" + firstdate
    outpng = plotFileDir + fileTag + ".png"
    plt.savefig(outpng,bbox_inches='tight')
    outpdf = plotFileDir + fileTag + ".pdf"
    plt.savefig(outpdf,bbox_inches='tight')
    print( "Wrote files %s and %s" % (outpng, outpdf))
else:
    # else show the plots
    plt.show()
