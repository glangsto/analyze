#Python Script to plot Calibrated NSF Horn Observations
#HISTORY
#19SEP11 GIL do not write the Kelvins file until fixed (.kel)
#19JUN29 GIL debugginginfo added
#18DEC11 GIL add argument processing loop, saving products
#18DEC10 GIL initial version based on m.py
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import radioastronomy
import hotcold
import copy
import interpolate

# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False # flag interpolate over spike in center of spectrum
doDebug = False   # flag printing debug info
#doDebug =True     # flag printing debug info
doSave = False    # flag saving intermediate files
flagRfi = True    # flag flagging RFI
doFold = False    # fold spectra to address an old issue; not normally used.
doSub = False     # define subtract baseline
doCalibrate = True# flagCalibrating observations
outDir = "./"     # define output directory for saving files
mytitle = ""      # define plot title
doHanning = False # hanning smooth hot load to reduce calibration noise

# put your list of known RFI features here.  Must have at least two, if flagRfi is true.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]   # integer number of channels to interpolate over
# currently used velocities for baseline fitting
maxvel = 300.
minvel = -300.
thot = 285.0  # define hot and cold in Kelvins
tcold = 10.0
nuh1 = 1420.40557177667 # neutral hydrogen frequency (MHz)
nureference = 1.E6*nuh1 # neutral hydrogen frequency (Hz)

nargs = len(sys.argv)
if nargs < 3:
    print "C: Calibrate Science Aficonado (NSF) horn observations"
    print "Usage: C [options]  <average_seconds> <files>"
    print ""
    print "Where many parameters are optional:"
    print "-B Subtract a linear baseline fit to Spectra at Min and Max Velocities"
    print "   Min and max default velocities are: %7.1f, %7.1f km/sec" % (minvel, maxvel)
    print "-C Flag the Center channel, using interpolation."
    print "   This removes a strong narrow feature created by many Software Defined Radios (SDRs)"
    print "-D Additional Debug printing."
    print "-I Optionally Flag Known Radio Frequency Interference (RFI)"
    print "   Note you need to update the c.py program to add your list of Frequencies"
    print "   RFI frequencies are Location Dependent"
    print "-H <hot load Temperature> Set the effective temperature of the hot load (Kelvins)"
    print "-N Not Calibrate.  This mode is used for tests of raw spectra"
    print "-O <output directory> Set the output directory for saved files"
    print "-R <Reference Frequency> Rest Frequency (Hz) used for Doppler calculations: %8.3f (MHz)" % (nureference*1.E-6)
    print "-S Save average spectra in files.  The Hot and Cold Load averages are saved, too."
    print "   Average spectra have -ave added to their names"
    print "   Calibrated spectra have a .kel (for Kelvins) extension"
    print "-X Hanning smooth the hot load observation to reduce calibration noise"
    print "-T <plot title String> Label for plot"
    print "-VA <low velocity> limit to exclude for baseline fitting"
    print "-VB <high velocity> limit to exclude for baseline fitting"
    print "Where:"
    print "   <average_seconds>: Number of seconds of observations to average."
    print "   <average_seconds> is clock time, not observing time, so 3600. gives one plot for each hour"
    print "   <files> are Horn Observation files"
    print "   <files> must include both data pointed up (.ast) and down (.hot) observations"
    print "      All .hot files are assumed to have a system temperature of %7.1f K" % thot
    print ""
    print " -- Glen Langston (glangsto@nsf.gov), 2018 December 11"
    print ""
    exit()

#first argument is the averaging time in seconds
timearg = 1
namearg = 2
iarg = 1          # start searching for input flags
# for all arguments, read list and exit when no flag argument found
while iarg < nargs:

    # if folding data
    if sys.argv[iarg].upper() == '-F':
        print 'Folding specectra'
        doFold = True
    elif sys.argv[iarg].upper() == '-I':
        flagRfi = True
    elif sys.argv[iarg].upper() == '-N':   # if no calibration 
        print 'Not calibrating Average Spectra in -ave files'
        doCalibrate = False
    elif sys.argv[iarg].upper() == '-C':
        flagCenter = True
    elif sys.argv[iarg].upper() == '-D':
        print 'Adding Debug Printing'
        doDebug = True
    elif sys.argv[iarg].upper() == '-S':
        print 'Saving Average Spectra in -ave files'
        doSave = True
    elif sys.argv[iarg].upper() == '-B':
        print 'Baseline subtraction'
        doSub = True
    elif sys.argv[iarg].upper() == '-X':
        print 'Hanning Smooth the Hot load Observation'
        doHanning = True
    elif sys.argv[iarg].upper() == '-O':   # now look for flags with arguments
        iarg = iarg+1
        timearg = timearg+1
        namearg = namearg+1
        outDir = sys.argv[iarg]
        print 'Output directory: ', outDir
    elif sys.argv[iarg].upper() == '-H':   # now look for flags with arguments
        iarg = iarg+1
        timearg = timearg+1
        namearg = namearg+1
        thot = float(sys.argv[iarg])
        print 'Hot Load Temperature: %7.2f K ', thot
    elif sys.argv[iarg].upper() == '-R':   # reference frequency (Hz)
        iarg = iarg+1
        timearg = timearg+1
        namearg = namearg+1
        nureference = float(sys.argv[iarg])
        print 'Reference Frequency: %8.3f (MHz): ', nureference * 1.E-6
    elif sys.argv[iarg].upper() == '-T':   # now look for flags with arguments
        iarg = iarg+1
        timearg = timearg+1
        namearg = namearg+1
        mytitle = sys.argv[iarg]
        print 'Plot title: ', mytitle
    elif sys.argv[iarg].upper() == '-VA':   # now look for flags with arguments
        iarg = iarg+1
        timearg = timearg+1
        namearg = namearg+1
        minvel = float(sys.argv[iarg])
        print 'Minimum velocity for baseline fit: %7.2f km/sec ' % (minvel)
        print '20 channels near this velocity will be used for the fit'
    elif sys.argv[iarg].upper() == '-VB':   # now look for flags with arguments
        iarg = iarg+1
        timearg = timearg+1
        namearg = namearg+1
        maxvel = float(sys.argv[iarg])
        print 'Maximum velocity for baseline fit: %7.2f km/sec ' % (maxvel)
        print '20 channels near this velocity will be used for the fit'
    else:
        break
    iarg = iarg + 1
    timearg = timearg+1
    namearg = namearg+1
# end of while not reading file names

# first argument is average time
try:
    avetimesec = float(sys.argv[timearg])
except:
    print 'First position argumen must be the average time in seconds'
    print 'First argument: ', sys.argv[timearg]
    exit()
avetime = datetime.timedelta(seconds=avetimesec)  # need datatime format

if doDebug:
    print "Average time: ", avetimesec, " (seconds)"

newObs = False
linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-c','-m','-y','-c','-m','-y','-c','-m','-y','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1.
xallmax = -9.e9
xallmin = 9.e9
yallmax = -9.e9
yallmin = 9.e9


tmin = 20.0  # Tsys never less than 20 K
tmax = 999.0 # define reasoanable value limits

velbaseline = [ minvel, maxvel] # array of velocities in km/sec

nplot = 0
nhot = 0         # number of obs with el < 0
ncold = 0
nhigh = 0        # highest galactic latitude
lastfreq = 0.
lastbw = 0.
lastgain = 0.
lastel = 0.
lastaz = 0.
firstdate = ""
lastdate = ""
firstaz = -1
otheraz = -1
dt = datetime.timedelta(seconds=0.)

#rest of arguments are the file names
names = sys.argv[namearg:]
names = sorted(names)
nFiles = len(names)

if doDebug:
    print 'First Name: ',names[0], namearg

hotnames = copy.deepcopy( names)
coldnames = copy.deepcopy( names)
plotnames = copy.deepcopy( names)

nhot, hot = hotcold.hotaverage( hotnames)
if doFold:
    hot.foldfrequency()

if nhot > 0:
    hot.ydataA = scalefactor * hot.ydataA 
    outstring = "Found %3d Hot Load observations          " % nhot
    sys.stdout.write(outstring)
    sys.stdout.flush()
else:
    print "No Hot Load observations, can not calibrate"
    exit()

nChan = hot.nChan
nData = nChan
n6 = int(nData/6)
n56= 5*n6
xv = hot.xdata[0:nChan] * 1.E-6
yv = copy.deepcopy( hot.ydataA[0:nChan])
if doDebug: 
    na = int(nChan/16)
    nb = int(15*nChan/16)
    print ''
    print 'Hot x,y[%4d]: %7.3f, %7.3f' % (na, xv[na], yv[na])
    print 'Hot x,y[%4d]: %7.3f, %7.3f' % (nb, xv[nb], yv[nb])
    print 'N Chan: ', nChan, len(xv)

# if flagging RFI
if flagRfi:
    hv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
else:
    hv = yv

if flagCenter:
    hv = hotcold.flagCenter( hv)

if doHanning:
    yv = copy.deepcopy(hv)
    for iii in range( 1,nData-1):
        hv[iii] = (yv[iii-1]+yv[iii]+yv[iii]+yv[iii+1])*.25

if doDebug:
    print "hv: %8.2f (counts)" % (hv[300])

if doSave:
    hot.writecount = 1
    hot.ydataA = hv
    hot.write_ascii_ave( outDir)
else:
    print ""
        
####################################### End of hot load processing
vel = hot.velocities(nureference)
chanbaseline = hot.vel2chan( velbaseline, nureference)
xa = int(chanbaseline[0])
xb = int(chanbaseline[1])

if xa > xb:   # keep the channel ranges from small to large
    temp = xa
    xa = xb
    xb = temp

if xa < 21:           # must avoid edge channels when fitting baselines 
    xa = 21
elif xa > nData - 21:
    xa = nData - 21

if xb > nData -  21:
    xb = nData - 21
elif xb < 21:
    xb = 21

if xa == xb:
    xa = xb - 1
xa = int(xa)
xb = int(xb)

if doDebug: 
    print ''
    print 'Hot File Coord: '
    print 'nChan, refChan: ',hot.nChan, hot.refChan
    print 'centerFreqHz  : ',hot.centerFreqHz
    print 'bandwidthHz   : ',hot.bandwidthHz
    print ''
    print 'N Channels    : ',nData
    print 'velocity range: ', velbaseline
    print 'Reference Freq. ', nureference
    print 'channel  range: ', chanbaseline
    print 'Min Vel at channel: ',xa, minvel
    print 'Max Vel at channel: ',xb, maxvel
                                   

ncold, cold, minel, maxel  = hotcold.coldaverage( coldnames)  # compute cold load average
if ncold < 1.:
    print "No high galactic latitude data: can not calibrate"
    exit()
else:
    cold.ydataA = scalefactor * cold.ydataA
    yv = copy.deepcopy(cold.ydataA)
    outstring = "Found %3d High Elevation spectra         " % (ncold)
    sys.stdout.write(outstring)
    sys.stdout.flush()
    
if doFold:
    cold.foldfrequency()

if flagRfi:
    cv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
else:
    cv = yv

if flagCenter:
    cv = hotcold.flagCenter( cv)


if doDebug:
    print "cv: %8.2f (counts)" % (cv[300])

if doSave:
    cold.writecount = 2
    cold.ydataA = copy.deepcopy(cv)
    cold.write_ascii_ave( outDir)
else:
    print ""

########################### End of Cold Spectrum
fig, ax1 = plt.subplots(figsize=(10, 6))

# compute gain in Kelvins/count
hotv = copy.deepcopy(hv)
Trx, gain = hotcold.compute_gain( hotv, cv, thot, tcold)

if doDebug:
    print "gv: %8.2f (Kelvins/count)" % (gain[300])
    print "Median Receiver + Antenna Temp: %8.2f (K)" % ( Trx)

if not doCalibrate: # if not calibrating, show hot and cold spectra
    tsky = copy.deepcopy(hv)
    if doSub:
        ya = np.median(tsky[(xa-10):(xa+10)])
        yb = np.median(tsky[(xb-10):(xb+10)])
        slope = (yb-ya)/(xb-xa)
        for iii in range( nData):
            tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))
    yallmin = min(tsky[xa:xb])
    yallmax = max(tsky[xa:xb])
    plt.plot(vel, tsky, '-r', linestyle='-', label="Hot  Load", lw=2)
    if doDebug:
        print 'hv[300]: ', tsky[300]

    tsky = copy.deepcopy(cv)
    if doSub:
        ya = np.median(tsky[(xa-10):(xa+10)])
        yb = np.median(tsky[(xb-10):(xb+10)])
        slope = (yb-ya)/(xb-xa)
        for iii in range( nData):
            tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))
    # reset the min/max values to include the hot and cold loads
    ymin = min(tsky[xa:xb])
    ymax = max(tsky[xa:xb])
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)

    plt.plot(vel, tsky, '-b', linestyle='-', label="Cold Load", lw=2)
    if doDebug:
        print 'cv[300]: ', tsky[300]

# end if not calibrating, then plot hot/cold

nRead = 0        
avenames = copy.deepcopy( plotnames)
# now read through all data and average cold sky obs
for filename in plotnames:
    parts = filename.split('/')
    nparts = len(parts)
    aname = parts[nparts-1]
    parts = aname.split('.')
    aname = parts[0]
    extension = parts[1]
    nRead = nRead + 1
# exclude hot load data for averaging
    if extension == 'hot':
        continue
# also exclude summaries
    if extension == 'sum':
        continue

    rs = radioastronomy.Spectrum()
#  print filename
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    if doFold:
        rs.foldfrequency()

    date, time = rs.datetime()

# if not a sky observation
    if rs.telel < 0.:
        continue
    if firstdate == "":            # if very first spectrum
        firstdate = date
        firstaz = rs.telaz
        lastaz = rs.telaz          # heep track of a range of azimuths
        otheraz = rs.telaz

    if otheraz != rs.telaz:
        otheraz = rs.telaz

    if minel > rs.telel:
        minel = rs.telel
    if maxel < rs.telel:
        maxel = rs.telel

# if first time reading data, set obs parameters
    if lastfreq == 0.:
        lastfreq = rs.centerFreqHz 
        lastbw = rs.bandwidthHz
        lastgain = rs.gains[0]
        lastaz = rs.telaz
        lastel = rs.telel
        cold = copy.deepcopy( rs)
        crosszero = False
        firstlon = rs.gallon
        ncold = 0

    if ncold > 0:
        # time difference is between mid-points of integrations
        dt = rs.utc - cold.utc 
        # add half of last observation
        dt = dt + datetime.timedelta(seconds=rs.durationSec/2.)
        # plus time before start of the first
        dt = dt + datetime.timedelta(seconds=cold.durationSec/2.)

        lastdate = date

        newAzEl = (lastaz != rs.telaz) or (lastel != rs.telel)
        newObs = (lastfreq != rs.centerFreqHz) or (lastbw != rs.bandwidthHz) or (lastgain != rs.gains[0]) or newAzEl
        if newObs:
            print "Change in observing parameters: "
            if lastfreq != rs.centerFreqHz:
                print "LastFreq: ", lastfreq/1e6, "New: ", rs.centerFreqHz/1e6, " MHz"
                lastfreq = rs.centerFreqHz
            if lastbw != rs.bandwidthHz:
                print "LastBandwidth: ", lastbw/1e6, "New: ", rs.bandwidthHz/1e6, " MHz"
                lastbw = rs.bandwidthHz
            if lastgain != rs.gains[0]:
                print "LastGain: ", lastgain, "New: ", rs.gains[0], " dB"
                lastgain = rs.gains[0]
            if newAzEl:
                print "LastAzEl: ", lastaz,lastel, "New: ", rs.telaz,rs.telel, " degrees"
                lastaz = rs.telaz
                lastel = rs.telel

        if nRead == nFiles:   # if this is the last file, must force output
            newObs = True

        # if time to average (or end of all files)
        if (dt > avetime) or newObs:
            nave, cold = hotcold.average( avenames[0:ncold])       # compute average of files
            if nave != ncold:
                print "Error averaging spectra", nave, ncold
            xv = copy.deepcopy(cold.xdata)
            xv = xv * 1.E-6                                  # convert to MHz
            yv = copy.deepcopy( cold.ydataA)                 # work with arrays
            yv = yv * scalefactor
            if flagRfi:
                yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi

            if flagCenter:
                yv = hotcold.flagCenter( yv)

            xmin = min(xv)
            xmax = max(xv)
            xallmin = min(xmin, xallmin)
            xallmax = max(xmax, xallmax)
            count = cold.count
            note = cold.noteA
            ncolor = min(nmax-1, nplot) 

            label = 'L,L=%5.1f,%5.1f' % (cold.gallon, cold.gallat)

            if doCalibrate: # if not calibrating, show hot and cold spectra
                tsky = hotcold.compute_tsky( xv, yv, gain, nureference)
                cold.bunit = 'Kelvins'
            else:
                tsky = yv
            cold.ydataA = copy.deepcopy(tsky)
            
            if doSub:
                ya = np.median(tsky[(xa-10):(xa+10)])
                yb = np.median(tsky[(xb-10):(xb+10)])
                slope = (yb-ya)/(xb-xa)
                for iii in range( nData):
                    tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))
                cold.ydataA = tsky

            ymed = np.median(cold.ydataA[n6:n56])
            ystd = np.std(cold.ydataA[n6:n56])
            if ystd <= 0.:
                ystd = 0.001
            ymin = min(cold.ydataA[xa:xb])
            ymax = max(cold.ydataA[xa:xb])

            if doDebug:
                print('Calibrated Max, Median and Std Dev: %8.3f %8.3f %8.3f' % (ymax, ymed, ystd))

            if doDebug:
                print 'xa,xb: ',xa,xb
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            date, time = cold.datetime()
            labeltime = date + " " + time
            if minel == maxel:
                label = '%s L,L=%5.1f,%4.1f' % (labeltime, cold.gallon, cold.gallat)
            else:
                label = '%s A,E=%.0f,%4.0f L,L=%5.1f,%4.1f' % (labeltime, cold.telaz, cold.telel, cold.gallon, cold.gallat)
            outstring = ' Max: %6.1f Med: %6.1f %c; %3d %s: ' % (ymax, ymed, cold.bunit[0], ncold, label)
            sys.stdout.write(outstring)
            sys.stdout.flush()

            if cold.gallat < 7.5 and cold.gallat > -7.5:
                plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
            elif cold.gallat < 15. and cold.gallat > -15.:
                plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
            else:
                plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
            nplot = nplot + 1
            ncold = 0
            if doSave:
                cold.writecount = nplot+2
                cold.write_ascii_ave( outDir)
            else:
                print ""
            if doDebug:
                print 'Vel, Tsky: %8.3f, %8.3f at %8.3f ' % (vel[300],tsky[300],xv[300])
                print 'Vel, Tsky: %8.3f, %8.3f at %8.3f ' % (vel[700],tsky[700],xv[700])
            # end if data to average and plot
        # end if a new observation

    if ncold <= 0:  # if starting a new observation
        cold = rs

    avenames[ncold] = filename
    ncold = ncold + 1
    # end if not a enough time

    # end if a cold file
    if nRead > nFiles:
        break

    #end for all files to sum
# end of all files to read

if doDebug:
    print 'Number of remaining observations not plotted: ', ncold

# if data remaing from summation
if ncold > 1:
    ncold, cold = hotcold.average( avenames[0:ncold-1])

    xv = copy.deepcopy(cold.xdata * 1.E-6)
    yv = cold.ydataA * scalefactor
    if flagRfi:
        yv = interpolate.lines( linelist, linewidth, xv, yv) # interpolate rfi
    
    if flagCenter:
        yv = hotcold.flagCenter( yv)
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin, xallmin)
    xallmax = max(xmax, xallmax)
    count = cold.count
    note = cold.noteA
    ncolor = min(nmax-1, nplot) 

    if doCalibrate: 
        # finally compute tsys
        tsky = hotcold.compute_tsky( xv, yv, gain, nureference)
        cold.bunit = 'Kelvins'
    else:
        tsky = yv

    if doSub:
        ya = np.median(tsky[(xa-10):(xa+10)])
        yb = np.median(tsky[(xb-10):(xb+10)])
        slope = (yb-ya)/(xb-xa)
        for iii in range( nData):
            tsky[iii] = tsky[iii] - (ya + (slope*(iii-xa)))
    cold.ydataA = copy.deepcopy( tsky)

    ymed = np.median(tsky)

    ymin = min(tsky[xa:xb])
    ymax = max(tsky[xa:xb])
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    date, time = cold.datetime()
    labeltime = date + " " + time
    if minel == maxel:
        label = '%s L,L=%5.1f,%4.1f' % (labeltime, cold.gallon, cold.gallat)
    else:
        label = '%s A,E=%.0f,%4.0f L,L=%5.1f,%4.1f' % (labeltime, cold.telaz, cold.telel, cold.gallon, cold.gallat)
    outstring = ' Max: %6.1f Med: %6.1f %c; %3d %s: ' % (ymax, ymed, cold.bunit[0], ncold, label)
    sys.stdout.write(outstring)
    sys.stdout.flush()
    if not doCalibrate:   # if not calibratiing just plot 
        tsky = yv
    if cold.gallat < 7.5 and cold.gallat > -7.5:
        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=4)
    elif cold.gallat < 15. and cold.gallat > -15.:
        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label, lw=2)
    else:
        plt.plot(vel, tsky, colors[ncolor], linestyle=linestyles[ncolor],label=label)
    nplot = nplot + 1

    if doSave:
        cold.writecount = nplot+2
        cold.write_ascii_ave( outDir)
    else:
        print ""
    if doDebug:
        print 'Vel, Tsky: ',vel[300],tsky[300],xv[300]
        print 'Vel, Tsky: ',vel[700],tsky[700],xv[700]
        
# if observations cover multiple days
if firstdate != lastdate:
    date = firstdate + " - " + lastdate
else:
    date = firstdate

if mytitle == "":
    mytitle = rs.site

if (firstaz == otheraz):
    mytitle = mytitle + "; Az = %6.1f, " % (firstaz)
else:
    mytitle = mytitle + "; Az = %6.1f to %6.1f, " % (firstaz, otheraz)

if minel == maxel:
    mytitle = mytitle + " El=%6.1f" % (minel)
else:
    mytitle = mytitle + " El=%6.1f to %6.1f" % (minel, maxel)

fig.canvas.set_window_title(mytitle)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_rotation('vertical')
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_rotation('vertical')
# select the relevant range of velocities (km/sec) for plotting
plt.xlim(minvel, maxvel)
# keep the plot from becoming too narrow
if yallmax < 8:
    yallmax = 8
# set the y scale to go above and below all data
plt.ylim((yallmin*.97)-1., 1.15*yallmax)
plt.title(mytitle, fontsize=16)
plt.xlabel('Velocity (km/sec)', fontsize=16)
if doCalibrate:
    plt.ylabel('Intensity (Kelvins)', fontsize=16)
else:
    plt.ylabel('Intensity (Counts)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
