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

class Plot(object):
    """
    Define the different types of plots
    """

#   some SDRs put spike in center of spectrum; indicate spike flagging here
    def __init__( self, verbose=True,  
                  flagCenter = False, 
                  doDebug = False, 
                  doSave = False, 
                  flagRfi = True, 
                  plotFrequency = True, 
                  myTitle = "",
                  doPlotFile = False, 
                  plotFileDir = "~/"  
                  ):
        self.verbose = verbose          # print extra messages
        self.flagCenter = flagCenter    # if flag spike in center of spectrum,
        self.doDebug = doDebug          # if printing debug info
        self.doSave = doSave            # if saving intermediate files
        self.flagRfi = flagRfi          # flag flagging RFI
        self.plotFrequency = plotFrequency   # plot frequency not velocity
        self.myTitle = myTitle          # provide plot title
        self.doPlotFile = doPlotFile    # if creating a plot file
        self.plotFileDir = plotFileDir  # location of optional plot file

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

        self.myTitle = myTitle

        # currently used velocities for plotting range
        self.maxvel = 180.
        self.minvel = -180.
        self.maxPlot = int(25)
        self.fileTag = ""
        self.firstdate = ""
        self.lastdate = ""
        self.doScaleAve = False
        self.xa = -1     # initialize indices corresponding to bandwidth to plot
        self.xb = -1
        self.scalefactor = 1.0
        self.nplot = 0
        return
    # end of init
    
    def help(self, argstring=""):
        ###
        # provide user help info also parses and updates configuration
        ###
        args = argstring.split()
        nargs = len(args)

        if (nargs < 1): 
            print("raw.help(flags): Help plotting raw counts of telescope observations")
            print("Usage: raw <flags> <files>")
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
                self.plotFileDir = args[iarg]
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
            elif args[iarg].upper() == '-CEN':   # if nU ref is provided (in MHz)n
                iarg = iarg+1
                self.nuRefFreq = float(args[iarg])
                print( 'Reference Frequency : %9.3f MHz' % (self.nuRefFreq))
            elif args[iarg].upper() == '-Z':     # label written files
                iarg = iarg+1
                self.fileTag = str(args[iarg])
                print( 'File tag: %s' % (self.fileTag))
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

    def raw(self, names):
        # to create plots in cronjobs, must use a different backend
        if self.doPlotFile:
            mpl.use('Agg')
        import matplotlib.pyplot as plt
        import gainfactor as gf

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

        names = names.split()
        nnames = len(names)

        print ("Raw: Plotting %d files" % (nnames))
        
        # plot no more than N spectra
        for filename in names:
            if len(filename) < 2:
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
    
            vel = np.zeros(nData)
            for jjj in range (0, nData):
                vel[jjj] = c * (self.nuRefFreq - xv[jjj])/self.nuRefFreq

            if not self.plotFrequency:
                self.xa, self.xb = gf.velocity_to_indicies( vel, minvel, maxvel)
    
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
        # if writing files
        if self.doPlotFile:
            self.firstdate = self.firstdate[2:]
            if self.fileTag == "":
                self.fileTag = "R-" + self.firstdate
            outpng = self.plotFileDir + self.fileTag + ".png"
            plt.savefig(outpng,bbox_inches='tight')
            outpdf = self.plotFileDir + self.fileTag + ".pdf"
            plt.savefig(outpdf,bbox_inches='tight')
            print( "Wrote files %s and %s" % (outpng, outpdf))
        else:
            # else show the plots
            plt.show()

        #end of ras.Raw(names)
        return

    def tsys(self, names):
        ###
        # Plot Tsys calibrated spectra
        ###
        print("Names:", (names))

        # End of ras.tsys()
        return
    
