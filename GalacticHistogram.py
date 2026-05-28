#/bin/python
#Calculate histogram of galactic latitude ranges, weighted by the time in range
#optionally plot histogram
#to test use:  python GalacticLatitudeHistogram.py data
#HISTORY
#26May28 Gil add option for event histograms in plot, remove ReadGalactic
#26May25 Gil annotate the plot with telescope and elevation
#26May24 Gil initial version
import sys
import numpy as np
import matplotlib.pyplot as plt
import fullpath

# types for histograms
NONTYPE = 0
ASTTYPE = 1
EVETYPE = 2
HOTTYPE = 3

class GalacticLatitudeHistogram:
    def __init__(self, bin_width_deg=5):
        self.bin_width = bin_width_deg
        self.edges = np.arange(-90, 90 + bin_width_deg, bin_width_deg)
        self.counts = np.zeros(len(self.edges) - 1, dtype=float)
        self.total_time_sec = 0.0

    def add_measurement(self, gal_lat_deg, duration_sec):
        """
        gal_lat_deg : float
            Galactic latitude in degrees.
        duration_sec : float
            Duration of this measurement in seconds.
        """
        idx = np.searchsorted(self.edges, gal_lat_deg, side="right") - 1
        if 0 <= idx < len(self.counts):
            self.counts[idx] += duration_sec
            self.total_time_sec += duration_sec

    def normalized(self):
        """
        Returns histogram normalized by total observing time,
        expressed as fraction of a day.
        """
        if self.total_time_sec == 0:
            return np.zeros_like(self.counts)

        seconds_per_day = 86400.0
        return self.counts / seconds_per_day

    def bin_centers(self):
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    def plot(self, title="Galactic Latitude Histogram", date="", show=True, \
             plotType=ASTTYPE):
        """
        plot() allows shows galactic distribution and
        diagnosis of issues with observations
        """
        centers = self.bin_centers()
        # convert fraction of day to minutes
        if plotType == ASTTYPE: 
            norm_hist = self.normalized() * 1440.
        else:
            norm_hist = self.normalized()
            
        plt.figure(figsize=(10, 6))
        plt.bar(centers, norm_hist, width=self.bin_width, \
                align="center", edgecolor="black")
        plt.xlabel("Galactic Latitude (deg)", fontsize=18)
        if plotType == EVETYPE:
            plt.ylabel("Event Signal/Noise Ratio",fontsize=18)
        else:
            plt.ylabel("Observing Time (Minutes)",fontsize=18)
        # special galactic latitude range for northern observations
        #        plt.xlim(-75., 90.)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title(title, fontsize=20)
        histmax = max( norm_hist)
        # annoate plot with date
        plt.text( -80., histmax*.95, date, fontsize=18)
        plt.grid(True, linestyle="--", alpha=0.5)

        # print("Check of Scaling: %7.1f" % (norm_hist.sum()))
        if show:
            plt.show()
    # end of plot
        return centers, norm_hist

def createLogName( logDirectory, date, telIndex):
    # assume date string has format 26May28
    # will separate measurement logs into years, months then days
    dateYear = str(date[0:2])    # two digits of year in date string: 26
    yearMonth = str(date[0:5])   # get Year+Month:  26May
    # create /home/karl/daily/2026/26May/26May28 string
    fullPath = logDirectory + "/20" + dateYear + "/" + yearMonth + "/" + date
    # create the directories leading to the log file
    fullpath.fullpath( fullPath)
    # finally prepare to write log file
    fullPath = fullPath + ("/spectraHistogram+%02d.log" % telIndex)
    return fullPath

def writeHistogram( logDirectory, date, dataDir, telIndex, degree, \
                    centers, norm_hist):
    """
    Write out the histogram of Observation Durations verse galactic Latitudes
    """
    fullPath = createLogName( logDirectory, date, telIndex)
    f = open( fullPath, 'w')
    f.write("#DATE     = %s \n" % (date))
    f.write("#TELINDEX = %d \n" % (telIndex))
    f.write("#DATADIR  = %s \n" % (dataDir))
    f.write("#BINWIDTH = %8.1f / Degrees \n" % (degree))
    nHist = len( norm_hist)
    f.write("#Center   Duration \n")
    f.write("# (deg)     (min)  \n")
    for iii in range( nHist):
        f.write("%8.1f %8.2f \n" % \
                (centers[iii], norm_hist[iii]))
# end of writeHistogram()
    return

def readHistogram( dataDir):
    """
    read the histogram of Observation Durations verse galactic Latitudes
    """
    # open and read all lines
    f = open( dataDir, 'r')
    # read
    lines = f.readlines()
    f.close()
    centers = []
    samples = []
    telIndex = 0
    date = ""
    datadir = ""
    degrees = 4.
    # for all lines in file
    for oneline in lines:
        # remove extra spaces
        aline = " ".join(oneline.split())
        aChar = aline[0]
        lineparts = aline.split(" ")
        nparts = len(lineparts)
        #ignore blank, partial lines
        if nparts < 1:
            continue
        print(lineparts)
        if lineparts[0] == "#DATE":
            date = lineparts[2]
        elif lineparts[0] == "#TELINDEX":
            telIndex = int(lineparts[2])
        elif lineparts[0] == "#DATADIR":
            dataDir = lineparts[2]
        elif lineparts[0] == "#BINWIDTH":
            degree = float( lineparts[2])
        elif aChar == "#":
            continue
        if nparts == 2:
            centers.append( float(lineparts[0]))
            samples.append( float(lineparts[1]))
    # end for all lines
    # now convert to arrays
    centers = np.array(centers)
    samples = np.array( samples)
    print("centers  : %s" % (type(centers)))
    print("samples  : %s" % (type(samples)))
    return date, dataDir, telIndex, degree, centers, samples
# end of readHistogram()
