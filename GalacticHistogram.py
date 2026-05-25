#/bin/python
#Calculate histogram of galactic latitude ranges, weighted by the time in range
#optionally plot histogram
#to test use:  python GalacticLatitudeHistogram.py data
#HISTORY
#26May25 Gil annotate the plot with telescope and elevation
#26May24 Gil initial version
import sys
import numpy as np
import matplotlib.pyplot as plt

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

    def plot(self, title="Galactic Latitude Histogram", show=True):
        """
        plot() allows diagnosis of issues with observations
        """
        centers = self.bin_centers()
        # convert fraction of day to minutes
        norm_hist = self.normalized() * 1440.

        plt.figure(figsize=(10, 6))
        plt.bar(centers, norm_hist, width=self.bin_width, \
                align="center", edgecolor="black")
        plt.xlabel("Galactic Latitude (deg)", fontsize=18)
        plt.ylabel("Observing Time (Minutes)",fontsize=18)
        # special galactic latitude range for northern observations
        #        plt.xlim(-75., 90.)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title(title, fontsize=20)
        plt.grid(True, linestyle="--", alpha=0.5)

        # print("Check of Scaling: %7.1f" % (norm_hist.sum()))
        if show:
            plt.show()

#Optionally run a diagnotic test
def main():
    """
    Main executable testing histogram
    """
    import radioastronomy
    import nameToIndex
    import astropy
    
    ifile = 1
    degrees = 4.
    
    nargs = len(sys.argv)
    if nargs < 2:
        print('GalacticHistogram: Count fraction of time in Galactic Latitude rranges')
        print('usage: GalaceticHistogram [-D degrees] dir1')
        print('where: ')
        print('       -D  <number> degree width of histogram bin')
        print('       -V  optionally print verbose debugging info')
        print('  dir1     Directory to compute histogram of observations')
        exit()

    dir1 = sys.argv[ifile]
    telIndex = nameToIndex.getIndex( dir1)

    # creast a list of files in the input directory
    from os import listdir
    from os.path import isfile, join
    files1 = [f for f in listdir(dir1) if isfile(join(dir1, f))]
    
    print("Degree width of histogram bins: %0.1f " % (degrees))

    print("%5d Files in Directory: %s" % (len(files1), dir1))
#    print files

    rs = radioastronomy.Spectrum()

    obs1s = files1
    nFile = 0
    nHot = 0
    minel = 99.
    maxel = -99.
    minaz = 365.
    maxaz = -365.
    # now for all files
    for filename in files1:
        parts = filename.split(".")
        nparts = len(parts)
        if nparts < 2:   # if not fooo.eve type file name
            continue
        if parts[nparts-1] == "ast":
            obs1s[nFile] = filename
            nFile = nFile + 1
        elif parts[nparts-1] == "hot":
            nHot = nHot + 1
        else:
            continue

    glh = GalacticLatitudeHistogram( bin_width_deg=degrees)
    # now get event signficances
    iii = 0
    total = 0.
    for filename in obs1s:
        fullname = join(dir1, filename)
        rs.read_spec_ast( fullname, headerOnly=True)
        rs.azel2radec()    # compute ra,dec from az,el
        tint = rs.durationSec
        lat = rs.gallat
        total = total + tint
        # now add to histogram 
        glh.add_measurement( lat, tint)
        if 50 * int(iii/50) == iii:
            print("File %5d: az,el = %5.1f,%6.1f => %5.1f,%6.1f" % \
                  (iii, rs.telaz, rs.telel, rs.gallon, rs.gallat))

        # keep track of location of observations
        if rs.telel > maxel:
            maxel = rs.telel
        if rs.telel < minel:
            minel = rs.telel
        if rs.telaz > maxaz:
            maxaz = rs.telaz
        if rs.telaz < minaz:
            minaz = rs.telaz
    
        iii = iii + 1

    title = "Tel %2d: Galactic Latitude Histogram" % telIndex
    if maxaz == minaz:
        title = title + (" Az: %.1f" % maxaz)
    else:
        title = title + (" Az: %.1f to %.1f" % (minaz, maxaz))
    if maxel == minel:
        title = title + (" El: %.1f" % maxel)
    else:
        title = title + (" El: %.1f to %.1f" % (minel, maxel))
    # now plot histogram
    glh.plot( title = title)

    print(title)
    print("Fraction of time Observing %.1f %s" % \
          (glh.total_time_sec/864., "%"))
    return

if __name__ == "__main__":
    main()



    
