#Python class to plot raw NSF spectra.
#HISTORY
#21Jun07 GIL rr -> ras, complete class
#21Jun06 GIL initial version 
#
import os
import matplotlib as mpl
import numpy as np
import radioastronomy
import interpolate
import gainfactor as gf

class Plot(object):
    """
    Define the different types of plots
    """

#   some SDRs put spike in center of spectrum; indicate spike flagging here
    def __init__( self, verbose=True,  
                  flagCenter = False, 
                  doSave = False, 
                  flagRfi = True, 
                  plotFrequency = True, 
                  myTitle = "",
                  doPlotFile = False, 
                  outFileDir = "../",
                  aveTimeSec=3600.,
                  telIndex = 0,
                  doBaseline = False,
                  writeTsys = False,
                  writeKelvin = False,
                  ):
        self.verbose = verbose          # print extra messages
        self.flagCenter = flagCenter    # if flag spike in center of spectrum,
        self.doSave = doSave            # if saving intermediate files
        self.flagRfi = flagRfi          # flag flagging RFI
        self.plotFrequency = plotFrequency   # plot frequency not velocity
        self.myTitle = myTitle          # provide plot title
        self.doPlotFile = doPlotFile    # if creating a plot file
        self.outFileDir = outFileDir  # location of optional plot file
        self.aveTimeSec = aveTimeSec    # location of optional plot file
        self.telIndex = telIndex        # add telescope number to plots
        self.doBaseline = doBaseline    # optionally subtract a baseline
        self.doZero = False             # optionally add a zero line
        self.writeTsys = writeTsys      # write log if integrated intenities
        self.writeKelvin = writeKelvin  # write calibrated intensities
        self.doKeep = False             # optionally write+cold raw files
        self.lowel = 10.                # minimum elevation for cold obs.
        
        # constants
        self.EPSILON = 1.E-10           # define a small number
        # put your list of known RFI features here.
        #Must have at least two, if flagRfi is true.
        self.linelist = [1400.00, 1420.0]  # RFI lines in MHz
        self.linewidth = [ 5, 5] # number of channels to flag

        # define reference frequency for velocities (MHz)
        self.nuh1 = 1420.40575 # neutral hydrogen frequency (MHz)
        self.nuoh1= 1612.231   # OH line
        self.nuoh2= 1665.402   # OH line
        self.nuoh3= 1667.359   # OH Line
        self.nuoh4= 1720.530   # OH Line

        # select the frequency for plotting velocities
        self.nuRefFreq = self.nuh1

        self.myTitle = myTitle    # input plot Title

        # set default names for files for hot and cold loads
        # null means compute 
        self.hotFileName = ""
        self.coldFileName = ""

        # currently used velocities for plotting range
        self.maxvel = 220.
        self.minvel = -220.
        self.maxSvel = 150.   # max velocity for numerical integration 
        self.minSvel = - self.maxSvel # min velocity for numerical integration 
        self.maxPlot = int(25)
        self.fileTag = ""
        self.firstdate = ""
        self.lastdate = ""
        self.doScaleAve = False
        self.xa = -1     # initialize indices corresponding to bandwidth to plot
        self.xb = -1
        self.scalefactor = 1.0
        self.nplot = 0
        # initialize of range for search
        self.minGlat = +90.
        self.maxGlat = -90.
        self.minel = 200.
        self.maxel = -200.

        # create places to store hot and cold files
        self.ave_hot  = radioastronomy.Spectrum()
        self.ave_cold = radioastronomy.Spectrum()
        self.nhot = 0
        self.ncold = 0
        self.nData = 128
        self.hv = np.zeros(self.nData)
        self.cv = np.zeros(self.nData)
        self.gain = np.zeros(self.nData)
        self.vel = np.zeros(self.nData)
        return
    # end of init

    def help(self, argstring=""):
        ###
        # provide user help info also parses and updates configuration
        ###
        args = argstring.split()
        nargs = len(args)

        if (nargs < 1): 
            print("ras.help(flags): Plotting Inputs for telescope obs.")
            print("Usage: .help('<flags> <files>')")
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
            print("-0 optionally plot zero intensity line(s)")
            print("-MINEL optionally set the lowest elevation allowed for calibration obs (default %7.1f)" % (self.lowel))            
            print("")
            print("Glen Langston - NSF June 7, 2021")
            return ["",""]

        iarg = 0
        # must be some arguments, parse them
        while iarg < nargs:

            # if folding data
            if args[iarg].upper() == '-F':
                print('Folding specectra')
                self.doFold = True
            elif args[iarg].upper() == '-A':
                self.toScalAve = True
            elif args[iarg].upper() == '-AVE':
                iarg = iarg + 1
                self.aveTimeSec = int( args[iarg])
                print('Averaging for %6.1f seconds' % (self.aveTimeSec))
            elif args[iarg].upper() == '-B':   # if setting beginning sample
                iarg = iarg + 1
                self.xa = int( args[iarg])
                print('Plotting starting at channel: %4d' % (self.xa))
            elif args[iarg].upper() == '-C':
                self.flagCenter = True
            elif args[iarg].upper() == '-E':   # if setting ending sample
                iarg = iarg + 1
                self.xb = int( args[iarg])
            elif args[iarg].upper() == '-H':
                iarg = iarg+1
                self.maxvel = np.float( args[iarg])
                print('Maximum (high) velocity for sum: %7.2f km/sec'\
                      % (self.maxvel))
            elif args[iarg].upper() == '-L':
                iarg = iarg+1
                self.minvel = np.float( args[iarg])
                print('Minium (low)  velocity for sum: %7.2f km/sec' \
                      % (self.minvel))
            elif args[iarg].upper() == '-N':   # if number of spectra to plot
                iarg = iarg+1
                self.maxPlot = int(args[iarg])
                if self.maxPlot < 1:
                    print("Not Plotting")
                else:
                    print("Plot will have a maximum of %d spectra" \
                          % (self.maxPlot))
            elif args[iarg].upper() == '-P':
                self.doPlotFile = True
                iarg = iarg+1
                self.outFileDir = args[iarg]
            elif args[iarg].upper() == '-Q':
                self.plotFrequency = True
            elif args[iarg].upper() == '-R':
                self.flagRfi = True
            elif args[iarg].upper() == '-T':   # if plot title provided
                iarg = iarg+1
                self.myTitle = args[iarg]
                print('Plot Title : ', self.myTitle)
            elif args[iarg].upper() == '-V':   # default is plotting Frequency
                self.plotFrequency = False              # plot velocity
            elif args[iarg].upper() == '-VA':   # now look for flags with arguments
                iarg = iarg+1
                self.minvel = float(args[iarg])
                print('Minimum velocity for baseline fit: %7.2f km/sec ' \
                      % (self.minvel))
            elif args[iarg].upper() == '-VB':   # now look for flags with arguments
                iarg = iarg+1
                self.maxvel = float(args[iarg])
                print('Maximum velocity for baseline fit: %7.2f km/sec ' \
                      % (self.maxvel))
            elif args[iarg].upper() == '-CEN':   # if nU ref is provided (in MHz)
                iarg = iarg+1
                self.nuRefFreq = float(args[iarg])
                print( 'Reference Frequency : %9.3f MHz' % (self.nuRefFreq))
            elif sys.argv[iarg].upper() == '-MINEL':  # if min elevation
                iarg = iarg + 1
                self.lowel = float( sys.argv[iarg])
                print(( "Using elevations > %7.2f (d) for Cold load obs." \
                        % (self.lowel)))
            elif sys.argv[iarg].upper() == '-W':
                self.writeKelvin = True
            elif args[iarg].upper() == '-Z':     # label written files
                iarg = iarg+1
                self.fileTag = str(args[iarg])
                print( 'File tag: %s' % (self.fileTag))
            elif sys.argv[iarg].upper() == '-0':
                doZero = True
                print('Plotting zero intensity lines')
            else:
                break
            iarg = iarg + 1
        # returns the file names/directories
        if iarg >= nargs:
            return ["",""]
        names = args[(iarg-1):]
        return names
    # end of help/parsing arguments

    def Help(self, argstring):
        ###
        # Alias for help()
        ###
        names = self.help(argstring)
        return names

    def average_spec( self, ave_spec, in_spec, nave, firstutc, lastutc):
        """
        Averages spectra/weight functions by observing duration
        input/output
        ave_spec   spectrum structure containing weighted average
        input:
        in_spec    spectrum to be added to the average
        in/out:
        nave       Count of number of spectra averaged so far.
                   nave = 0 implies initializing the sum.
        firstutc   date of first observation averaged
        lastutc    date of current last date of spectrum to be averaged
        """

        if self.verbose:
            nData = len(in_spec.ydataA)
            n6 = int(nData/6)
            n56 = int(5*n6)            
            medianData = np.median( in_spec.ydataA[n6:n56])
            medianScale = np.median( in_spec.ydataA[n6:n56])
            print(( "Input: %8.3f, count: %6d; Scale: %8.3f" \
                    % (medianData, in_spec.count, medianScale)))

        # if restarting the sum
        if nave == 0:
            ave_spec = copy.deepcopy(in_spec)  # initial spectrum is 1st
            firstutc = in_spec.utc
            lastutc = in_spec.utc
            self.nData = len(in_spec.ydataA)
            nave = 1
            # replace with duration scaled intensities
            ave_spec.ydataA = (in_spec.ydataA * in_spec.durationSec)  
            # keep track of observing time for weighted sum
            ave_spec.durationSec = in_spec.durationSec
        else: # else not enough time yet, average ave_spec data
            if in_spec.utc < firstutc:
                firstutc = in_spec.utc
            elif in_spec.utc > lastutc:
                lastutc = in_spec.utc
            nave = nave + 1
            ave_spec.ydataA = ave_spec.ydataA + \
                (in_spec.ydataA * in_spec.durationSec)
            # keep track of observing time for weighted sum
            ave_spec.durationSec = ave_spec.durationSec + in_spec.durationSec
        
        return ave_spec, nave, firstutc, lastutc
    # END OF average_spec

    def normalize_spec( self, ave_spec, firstutc, lastutc):
        """
        Normaize the average after completing sum
        input/output
        ave_spec   input raw sum of observation 
               output normalized sum of observations
        input      first and last utcs in observation
        Note:  The count of observations is not used, rather the integration
        is weighted by durations.  This corrects for observations with 
        different integration times.
        """
        # now renormalize for total integration time
        ave_spec.ydataA = ave_spec.ydataA/float(ave_spec.durationSec)
        # compute average time from first and last utcs
        aveutc, duration = radioastronomy.aveutcs( firstutc, lastutc)
        ave_spec.utc = aveutc
        #need to re-calculate representative RA,Dec for average time
        ave_spec.azel2radec()

        # end of normalize_spec()
        return ave_spec

    def read_hot( self, names, ave_hot):
        """
        read_hot() reads in all files in the names list and averages hot load
        observations.   The hot load file is returned.
        While reading the files, the minmum elevation and 
        galactic latitudes are recorded
        """

        rs = radioastronomy.Spectrum()
        nhot = 0       # init count of hot files
        ncold = 0
        minGlat = +90.
        maxGlat = -90.
        minel = 200.
        maxel = -200.

        # if reading a hot file
        if (self.hotFileName != ""):
            self.ave_hot.read_spec_ast(self.hotFileName)
            if doScaleAve:
                ave_hot.ydataA = ave_hot.ydataA/ave_hot.count
            else:
                ave_hot.ydataA = ave_hot.ydataA/ave_hot.nave
            self.nhot = 1
            return
        
        # now run through and find hot and cold loads obs
        for filename in names:

            parts = filename.split('/')
            nparts = len(parts)
            if nparts == 1:
                aname = parts[0]
            else:
                aname = parts[nparts-1]
            parts = aname.split('.')
            nparts = len(parts)
            if nparts < 2:
                print(( 'File is not an astronomy file: %s' % (filename)))
                continue
            else:
                extension = parts[nparts-1]
            extension = extension.upper()
            if (extension != 'HOT') and (extension != 'AST') and \
               (extension != 'CLD'):
                print(( 'Extension not recognized : %s' % ( parts[nparts-1])))
                continue

            rs.read_spec_ast(filename)
            if self.doScaleAve:
                rs.ydataA = rs.ydataA/rs.count
            else:
                rs.ydataA = rs.ydataA/rs.nave
        
            nChan = len( rs.ydataA)
            if nChan != 32 and nChan != 64 and nChan != 128 and nChan != 256 \
               and nChan != 512 and nChan != 1024 and nChan != 2048 and \
                   nChan != 4096:
                print("Unusual data length (%d), file %s" % (nChan, filename))
                continue
            rs.azel2radec()    # compute ra,dec from az,el

            # if a hot load observation
            if rs.telel < 0:
                if nhot == 0:
                    firstutc = rs.utc
                    lastutc = rs.utc
                    minel = rs.telel

                # accumulate spectra, nhot is updated also
                self.ave_hot, nhot, firstutc, lastutc = \
                    average_spec( self.ave_hot, rs, nhot, firstutc, lastutc)
                if minel > rs.telel:
                    minel = rs.telel
            else: # else above horizon, find min, max galactic latitudes
                if ncold == 0:
                    minGlat = rs.gallat
                    maxGlat = rs.gallat
                    minGlon = rs.gallon
                    maxGlon = rs.gallon
                else:
                    if rs.gallat > maxGlat:
                        maxGlat = rs.gallat
                    if rs.gallat < minGlat:
                        minGlat = rs.gallat
                    if minel > rs.telel:
                        minel = rs.telel
                    if maxel < rs.telel:
                        maxel = rs.telel
                    ncold = ncold + 1
        # end for all files

        if nhot > 0:
            print(( "Found %3d Hot load observations" % (nhot)))
            self.ave_hot = normalize_spec( self.ave_hot, firstutc, lastutc)
        else:
            print( "No Hot load data, can not calibrate")
            return
            # exit()

        self.nhot = nhot
            
        # do more cleanup on spectra for RFI
        yv = copy.deepcopy(self.ave_hot.ydataA)
        if self.flagRfi:
            xv = ave_hot.xdata * 1.E-6
            # interpolate rfi
            hv = interpolate.lines( linelist, linewidth, xv, yv)
        else:
            hv = yv

        if self.flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
            icenter = int(self.nData/2)
            hv[icenter] = (hv[icenter-2] + hv[icenter+2])*.5
            hv[icenter-1] = (3.*hv[icenter-2] + hv[icenter+2])*.25
            hv[icenter+1] = (hv[icenter-2] + 3.*hv[icenter+2])*.25

        # all outputs are part of object
        self.ave_hot.ydataA = hv
        self.hv = hv
        self.minel = minel
        self.maxel = maxel
        self.minGlat = minGlat
        self.maxGlat = maxGlat
        
        # if keeping hot and cold files
        if self.doKeep and self.hotFileName == "":
            outname = radioastronomy.utcToName( ave_hot.utc)
            outname = outname + ".hot"  # output in counts
            # add telescope index
            outname = ("T%d-" % cpuIndex)+outname
            n = ave_hot.nChan
            print("Ave Hot: %d: %.6f" % \
                  (n, np.median( ave_hot.ydataA[int(n/3):int(2*n/3)])))
            ave_hot.write_ascii_file(self.outputDir, outname)
            print( "Wrote Average Hot  Load File: %s%s" % ("../", outname))

        if verbose:
            print("Min, Max Galactic Latitude: %7.1f,%7.1f" \
                   % (minGlat, maxGlat))
            print("Min, Max Elevation:         %7.1f,%7.1f" \
                   % (minel, maxel))

        # end of read hot
        return

    def read_angles( self, names):
        """
        read_angles() reads all files and counts number of files with 
        high elevation and high galactic latitude
        Inputs:
        lowel   minimum elevation to accept for cold load 
        """

        names = self.splitNames(names)
        
        # count types
        ncold = 0
        rs = radioastronomy.Spectrum()
        nName = len(names)
        minel = 90.
        maxel = -90.
        minGlat = 90
        maxGlat = -90

        # now average coldest data for calibration
        for filename in names:

            rs.read_spec_ast(filename)
            rs.azel2radec()    # compute ra,dec from az,el

            if rs.telel < self.lowel:  #if elevation too low for a cold load obs
                continue
            if rs.gallat > maxGlat:
                maxGlat = rs.gallat
            if rs.gallat < minGlat:
                minGlat = rs.gallat
            if minel > rs.telel:
                minel = rs.telel
            if maxel < rs.telel:
                maxel = rs.telel
        
            # count number of high elevation files
            if rs.telel > self.lowel:
                ncold = ncold + 1

        print(( "Found %d high elevation obs in %d files" % (ncold, nName)))
        self.ncold = ncold
        self.minel = minel
        self.maxel = maxel
        self.minGlat = minGlat
        self.maxGlat = maxGlat

        # end of read_angles()
        return ncold

    def read_cold( names, ave_cold, lowel, lowGlat):
        """
        read_cold() checks files and averages selected files with 
        high elevation and galactic Latitude.
        Inputs:
        lowel   minimum elevation to accept for cold load 
        lowGlat minimum galactic latitude
        """
        # starting new sum of cold (high elevation and galactic latitude) obs
        ncold = 0
        rs = radioastronomy.Spectrum()
        nName = len(names)

        # now average coldest data for calibration
        for filename in names:

            rs.read_spec_ast(filename)
            if rs.telel < self.lowel:  #if elevation too low for a cold load obs
                continue

            if doScaleAve:
                rs.ydataA = rs.ydataA/rs.count
            else:
                rs.ydataA = rs.ydataA/rs.nave
                rs.azel2radec()    # compute ra,dec from az,el
                
            # note this test excludes low galactic latitude ranges
            if rs.gallat > lowGlat or rs.gallat < -lowGlat:
                if ncold == 0:
                    firstutc = rs.utc
                    lastutc = rs.utc
                ave_cold, ncold, firstutc, lastutc = \
                    average_spec( ave_cold, rs, ncold, firstutc, lastutc)

            # end of all files to average
        if ncold < 1:
            print( "No high elevation data: can not calibrate")
            return 0
        else:
            ave_cold = normalize_spec( ave_cold, firstutc, lastutc)

        if flagRfi:
            cv = ave_cold.ydataA
            xv = ave_cold.xdata*1.E-6
            # interpolate rfi
            yv = interpolate.lines( linelist, linewidth, xv, cv)
        else:
            yv = ave_cold.ydataA
                
        if self.flagCenter:             # if flagging spike in center of plot
            # remove spike in center of the plot
            icenter = int(nData/2)
            yv[icenter] = (yv[icenter-2] + yv[icenter+2])*.5
            yv[icenter-1] = (3.*yv[icenter-2] + yv[icenter+2])*.25
            yv[icenter+1] = (yv[icenter-2] + 3.*yv[icenter+2])*.25
            
        ave_cold.ydataA = cv
        self.ave_cold = copy.deepcopy( ave_cold)
        self.cv = cv
        self.ncold = ncold
            
        if verbose:
            print( "Found %3d High Galactic Latitude spectra" % (ncold))
        return ncold

    def computeGain( self):
        """
        Compute the telescope gain (Kelvins per Count) based on hot and cold
        """
        
    def splitNames( self, names):
        """
        splitNames parses the names in a list and also expands directories
        to include all files in any directories found
        """
        names = names.split()
        names = sorted(names)
        nnames = len(names)
        print ("Found %d file names" % (nnames))

        # end of splitNames
        return names
        
    def raw(self, names):
        """
        Plot all files in the names list
        """
        # to create plots in cronjobs, must use a different backend
        if self.doPlotFile:
            mpl.use('Agg')
        import matplotlib.pyplot as plt

        if self.plotFrequency:
            print( "Ploting Intensity versus Frequency")
        else:
            print( "Ploting Intensity versus Velocity")

        linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
        colors = ['g', 'b', 'r', '-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

        xallmax = -9.e9
        xallmin =  9.e9
        yallmax = -9.e9
        yallmin =  9.e9
        c = 299792.458  # (Speed of light  km/sec)

        #for symbol, value in locals().items():

        # initialize spectrum for reading and plotting
        rs = radioastronomy.Spectrum()

        names = self.splitNames(names)
        
        # plot no more than N spectra
        for filename in names:
            # a file name must have a 3 letter extension (ie .hot, .ast)
            if len(filename) < 5:
                continue
            if self.verbose:
                print('%s' % (filename))

            rs.read_spec_ast( filename)
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

            if self.xa < 0:
                self.xa = 0
            if self.xb < 0:
                self.xb = nData
            if self.xb > nData:
                self.xb = nData
        
            if self.firstdate == "":
                self.firstdate = date
            self.lastdate = date
    
            self.vel = np.zeros(nData)
            for jjj in range (0, nData):
                self.vel[jjj] = c * (self.nuRefFreq - xv[jjj])/self.nuRefFreq

            if not self.plotFrequency:
                self.xa, self.xb = gf.velocity_to_indicies( self.vel, \
                                                            minvel, maxvel)
    
            # normize for different integration times
            if self.doScaleAve:
                rs.ydataA = rs.ydataA/rs.count
            else:
                rs.ydataA = rs.ydataA/rs.nave

            # must first determine if data are already scaled
            yv = rs.ydataA * self.scalefactor
            ymedian = np.median(yv[n6:n56])

            # if the latest software, the scale factor is just 1.
            if ymedian < .001:
                yv = rs.ydataA * self.scalefactor
            else:
                self.scalefactor = 1.0

            if not self.plotFrequency:
                xv = vel
            xmin = min(xv[self.xa:self.xb])
            xmax = max(xv[self.xa:self.xb])
            xallmin = min(xmin,xallmin)
            xallmax = max(xmax,xallmax)

            if self.flagRfi:
                # interpolate rfi
                yv = interpolate.lines( self.linelist, self.linewidth, xv, yv) 

            if self.flagCenter:             # if flagging spike in center of plot
                # remove spike in center of the plot
                icenter = int(nData/2)
                yv[icenter] = (yv[icenter-2] + yv[icenter+2])*.5
                yv[icenter-1] = (3.*yv[icenter-2] + yv[icenter+2])*.25
                yv[icenter+1] = (yv[icenter-2] + 3.*yv[icenter+2])*.25
            
            ymin = min(yv)
            ymax = max(yv)
            ymed = np.median(yv)

            print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' \
              % (ymax, ymed, ymax/ymed, rs.count, label))
            if self.nplot <= 0 and self.maxPlot > 0:
                fig,ax1 = plt.subplots(figsize=(10,6))
                fig.canvas.set_window_title(date)
                for tick in ax1.xaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 
                for tick in ax1.yaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 

            self.nplot = self.nplot + 1
            if self.nplot > self.maxPlot:
                break
            note = rs.site
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            plt.xlim(xallmin,xallmax)

            # scale min and max intensities for nice plotting
            plt.ylim(0.9*yallmin,1.25*yallmax)

            if self.plotFrequency:
                plt.plot(xv[self.xa:self.xb], yv[self.xa:self.xb], \
                         colors[self.nplot], \
                         linestyle=linestyles[self.nplot-1],label=label, lw=2)
            else:
                plt.plot(xv[self.xa:self.xb], yv[self.xa:self.xb], \
                         colors[self.nplot], \
                         linestyle=linestyles[self.nplot-1],label=label, lw=2)
        # end for all names
        if (self.maxPlot < 1) or (self.nplot < 1):
            print("No Plots, exiting")
            exit()
    
        if self.myTitle == "":
            self.myTitle = note
    
        plt.title(self.myTitle, fontsize=16)
        plt.xlabel('Frequency (MHz)',fontsize=16)
        ylabel = 'Intensity (%s)' % rs.bunit
        plt.ylabel(ylabel, fontsize=16)
        plt.legend(loc='upper right')
        # zero the plot count for next execution
        self.nplot = 0 
        # if writing files
        if self.doPlotFile:
            self.firstdate = self.firstdate[2:]
            if self.fileTag == "":
                self.fileTag = "R-" + self.firstdate
            outpng = self.outFileDir + self.fileTag + ".png"
            plt.savefig(outpng,bbox_inches='tight')
            outpdf = self.outFileDir + self.fileTag + ".pdf"
            plt.savefig(outpdf,bbox_inches='tight')
            print( "Wrote files %s and %s" % (outpng, outpdf))
        else:
            # else show the plots
            plt.show()

        #end of ras.raw(names)
        return

    def tsys(self, names):
        ###
        # Plot Tsys calibrated spectra
        ###
        names = self.splitNames(names)

        print("Names:")
        nnames = len(names)
        for iii in range(min(nnames,4)):
            print("%4d: %s" % (iii, names[iii]))

        # zero the plot count for next execution
        self.nplot = 0 

        # End of ras.tsys()
        return
    
