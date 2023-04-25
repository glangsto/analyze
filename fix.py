#Python Script to fix set of NSF events or spectra
#The list of header items to fix is given in help.
#These include El (elevation) and Az (azimuth)
#HISTORY
#23APR25 GIL if fixing center frequency also update xdata
#23MAR31 GIL add count of files fixed, add gain1
#21APR09 GIL deal with note files add fix telescope altitude
#19NOV22 GIL don't plot spectra
#19SEP24 GIL turn off plotting by default
#19FEB20 GIL initial version
#
import matplotlib.pyplot as plt
import sys
import radioastronomy
import interpolate
import numpy as np
import os

dy = -1.

nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

scalefactor = 1e8
xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9
# set magic values to identify new inputs
NOVALUE = -200.
newEl = NOVALUE
newAz = NOVALUE
newdEl = NOVALUE
newdAz = NOVALUE
newLat = NOVALUE
newLon = NOVALUE
newAlt = NOVALUE
newCen = NOVALUE
newBan = NOVALUE
newRefSample = NOVALUE
newRefChan = NOVALUE
newNChan = NOVALUE
newNTime = -200
newGain1 = NOVALUE
newGain2 = NOVALUE
newGain3 = NOVALUE
observer = ""
note = ""
telescope = ""
device = ""
# flag replacing file
replace = False
aFix = False
doPlot = False

ifile = 1
iii = ifile

while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-EL":
        newEl = float( sys.argv[iii+1])
        iii = iii + 1
        print("New El: %7.2f" % (newEl))
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:4]) == "-DEL":
        newdEl = float( sys.argv[iii+1])
        iii = iii + 1
        print("New El Offset: %7.2f" % (newdEl))
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-AZ":
        newAz = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Az: %7.2f" % (newAz))
        ifile = ifile + 2
        aFix = True
    if str(anarg[0:4]) == "-DAZ":
        newdAz = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Az Offset: %7.2f" % (newdAz))
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-LA":
        newLat = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Latitude : %13.9f" % (newLat))
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-LO":
        newLon = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Longitude: %13.9f" % (newLon))
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-AL":
        newAlt = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Altitude : %6.3f" % (newAlt))
        ifile = ifile + 2
        aFix = True
    if anarg == "-NT":
        newNTime = np.int( sys.argv[iii+1])
        iii = iii + 1
        print("New Number of Samples: ", newNTime)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-OB":
        observer = sys.argv[iii+1]
        iii = iii + 1
        print("Observers: ", observer)
        ifile = ifile + 2
        aFix = True
    if anarg[0:2] == "-P":
        doPlot = True
        iii = iii + 1
        print("Plotting fixed observations")
        ifile = ifile + 1
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print("Note: ", note)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-TE":
        telescope = sys.argv[iii+1]
        iii = iii + 1
        print("Telescope: ", telescope)
        ifile = ifile + 2
        aFix = True
    if anarg[0:4] == "-DEV":
        device = sys.argv[iii+1]
        iii = iii + 1
        print("Device: ", device)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-CE":
        newCen = float(sys.argv[iii+1])
        iii = iii + 1
        print("Center Frequency: ", newCen, " (MHz)")
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-BA":
        newBan = float(sys.argv[iii+1])
        iii = iii + 1
        print("Bandwidth ", newBan, " (MHz)")
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-RC":
        newRefChan = float(sys.argv[iii+1])
        iii = iii + 1
        print("Ref Chan  ", newRefChan)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-NC":
        newNChan = float(sys.argv[iii+1])
        iii = iii + 1
        print("N Channel ", newNChan)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-RS":
        newRefSample = float(sys.argv[iii+1])
        iii = iii + 1
        print("Reference Sample ", newRefSample, " ")
        ifile = ifile + 2
        aFix = True
    if anarg[0:6] == "-GAIN1":
        newGain1 = float(sys.argv[iii+1])
        iii = iii + 1
        print("Gain 1: %7.2f" % (newGain1))
        ifile = ifile + 2
        aFix = True
    if anarg[0:6] == "-GAIN2":
        newGain2 = float(sys.argv[iii+1])
        iii = iii + 1
        print("Gain 2: %7.2f" % (newGain2))
        ifile = ifile + 2
        aFix = True
    if anarg[0:6] == "-GAIN3":
        newGain3 = float(sys.argv[iii+1])
        iii = iii + 1
        print("Gain 3: %7.2f" % (newGain3))
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-RE":
        replace = True 
        print("Replacing original file")
        ifile = ifile + 1
        aFix = True
    iii = iii + 1

#print "Afix: ",aFix, iii
# if nothing to fix, give help
if aFix == False:
    print("FIX: Fix observing file parameters")
    print("Usage: Fix [-el elevation] [-az azimuth]... <file 1> [<file 2>] ... [<file N>]")
    print("Where optionally the following paramters may be fixed")
    print(" -az  Telescope azimuth in degrees")
    print(" -el  Telescope elevation in degrees")
    print(" -daz Telescope azimuth offset in degrees")
    print(" -del Telescope elevation offset in degrees")
    print(" -lat Telescope latitude in degrees")
    print(" -lon Telescope longitude in degrees")
    print(" -alt Telescope altitude above sea level in meters")
    print(" -cen Center Frequency (MHz)")
    print(" -ban Bandwidth (MHz)")
    print(" -rch Reference channel  (1 based)")
    print(" -nch Number of channels (1 based)")
    print(" -obs Observers names (ascii)")
    print(" -tel Telescope name (ascii)")
    print(" -not Note describing observation (ascii)")
    print(" -dev Software Defined Radio used for observation (ascii)")
    print(" -nt  Number of time samples in the observations")
    print(" -rs  Index of the Reference Time Sample")
    print(" -re  Replace original file with revised header")
    print(" -p   Plot fixed events")
    exit()

nplot = 0
nFix = 0
nfiles = nargs-ifile
for iii in range(nfiles):

    filename = sys.argv[iii+ifile]

    rs = radioastronomy.Spectrum()
#    print filename
    rs.read_spec_ast( filename)
    if newAz != NOVALUE:
        rs.telaz = newAz
    if newEl != NOVALUE:
        rs.telel = newEl
    if newdAz != NOVALUE:
        rs.teldaz = newdAz
    if newdEl != NOVALUE:
        rs.teldel = newdEl
    if newLat != NOVALUE:
        rs.tellat = newLat
    if newLon != NOVALUE:
        rs.tellon = newLon
    if newAlt != NOVALUE:
        rs.telelev = newAlt
    if newGain1 != NOVALUE:
        rs.gains[0] = newGain1
    if newGain2 != NOVALUE:
        rs.gains[1] = newGain2
    if newGain3 != NOVALUE:
        rs.gains[2] = newGain3
    if newAlt != NOVALUE:
        rs.telelev = newAlt
    if observer != "":
        rs.observer = observer
    if device != "":
        rs.device = device
    if note != "":
        rs.noteA = note
    if telescope != "":
        rs.site = telescope
    if newNTime > 0:
        rs.nTime = newNTime
    if newBan > 0:
        rs.bandwidthHz = newBan * 1.E6
    if newCen > 0:
        rs.centerFreqHz = newCen * 1.E6
    if newRefSample != NOVALUE:
        rs.refSample = newRefSample
    if newRefChan != NOVALUE:
        rs.refChan = newRefChan
    if newNChan != NOVALUE:
        rs.nChan = newNChan

    rs.azel2radec()    # compute ra,dec from az,el and telescope location

    # if changing center frequency or bandwidth, must also update X data
    if newCen > 0 or newBan > 0:
        nData = len(rs.xdata)
        dX = rs.bandwidthHz / float(nData)
        X0 = rs.centerFreqHz - (rs.refChan*dX)
        for iii in range( rs.nChan):
            rs.xdata[iii] = X0
            X0 = X0 + dX
    # now prepare to write
    parts = filename.split('/')
    nparts = len(parts)
    # get the file name without directory name
    aname = parts[nparts-1]
    filepart = aname
    # if no directory name
    if nparts == 1:
        dirname = "./"   # use current directory
    else:
        dirname = ""
        for i in range(nparts-1):
            dirname = dirname + parts[i] + "/"
#    print "Directory: ", dirname

    strutc = rs.utc.isoformat()
    parts = strutc.split('.')
    aname = parts[0]
    parts = aname.split('T')
    # the date and time are only used for labeling plots
    date  = parts[0]
    time  = parts[1]
    time  = time.replace('_',':')

    extension = ""
    # if replacing original file
    if replace:
        try:
            os.remove(filename)
        except:
            print("Cound not remove file: ",filename)
        outname = filename
        parts = filename.split('.')
        # last bit of file name is the extension
        nparts = len(parts)
        extension = parts[nparts-1]
    else:
        parts = filepart.split('.')
        nparts = len(parts)
        if (nparts == 2):
            extension = parts[1]
            outname = parts[0] + "-fix." + extension
        elif (nparts == 1):
            outname = parts[0] + "-fix"
        else:
            extension = parts[nparts-1]
            outname = parts[0] + "-fix." + extension
    
# now if a spectrum, fix name for elevation above zero 
# Spectra do not have time series.
    if rs.nTime <= 0:
        parts = outname.split('.')
        nparts = len(parts)
        # last part of file name is file type
        filetype = parts[nparts-1]
        # special case of a notes files
        if extension == "not":
            filetype = extension
        else:
            if rs.telel > 0.:
                filetype = 'ast'
            else:
                filetype = 'hot'
                
        if nparts == 2:
            outname = parts[0] + "." + filetype
        elif nparts == 1:
            outname = parts[0]
        else:  # else multiple name parts separated by "." 
            outname = ""
            for iii in range(nparts-1):
                outname = outname + parts[iii] + "."
            outname = outname + filetype
#    print "Output file name: ", outname
    rs.write_ascii_file( dirname, outname)
    nFix = nFix + 1

# only plot the first few events
    gallon = rs.gallon
    gallat = rs.gallat
    label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % ( time,rs.telaz,rs.telel,gallon,gallat)
    xs = rs.xdata 
    ya = rs.ydataA
    yb = rs.ydataB

# plotting only works for events
    if rs.nTime < 1:
        continue

    if nplot > 10:
        continue
    if not doPlot:
        continue
    nplot = nplot+1

    # the remainder of this code only plots events.
    xv = np.zeros(rs.nSamples*2)
    yv = np.zeros(rs.nSamples*2)
    j = 0
    dt = 0.5/rs.bandwidthHz
    t = xs[0]
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
    ymed = (ymin+ymax)/2.
    count = rs.count

    if nplot <= 0:
        fig,ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(date)
        nplot = nplot + 1

    print((' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label)))
    note = rs.noteA
#    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    plt.xlim(xallmin,xallmax)
#    plt.ylim(0.9*ymin,1.5*yallmax)
#    plt.ylim(0.9*ymin,1.25*yallmax)
    plt.ylim(yallmin,1.25*yallmax)

    plt.plot(xv, yv, colors[nplot], linestyle=linestyles[nplot],label=label)
# put labels on plot
if doPlot:
    plt.title(note)
    plt.xlabel('Time (s)')
    plt.ylabel('Intensity (Counts)')
    plt.legend(loc='upper right')
    plt.show()

print("Fixed %d Files" % (nFix))
