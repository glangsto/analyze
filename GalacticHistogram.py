#/bin/python
#Calculate histogram of galactic latitude ranges, weighted by the time in range
#HISTORY
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
        centers = self.bin_centers()
        norm_hist = self.normalized()

        plt.figure(figsize=(10, 5))
        plt.bar(centers, norm_hist, width=self.bin_width, \
                align="center", edgecolor="black")
        plt.xlabel("Galactic Latitude (deg)")
        plt.ylabel("Normalized Time (fraction of day)")
        plt.title(title)
        plt.grid(True, linestyle="--", alpha=0.5)

        if show:
            plt.show()

def main():
    """
    Main executable for gridding astronomical data
    """
    import radioastronomy
    import astropy
    
    ifile = 1
    degrees = 5.
    
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
    print("Dir 1: ", dir1)

    from os import listdir
    from os.path import isfile, join
    files1 = [f for f in listdir(dir1) if isfile(join(dir1, f))]
    
    print("Degree width of histogram bins: %0.1f " % (degrees))

    print("%5d Files in Directory 1" % (len(files1)))
#    print files

    rs = radioastronomy.Spectrum()

    obs1s = files1
    nFile = 0
    nHot = 0

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
        rs.read_spec_ast( fullname)
        rs.azel2radec()    # compute ra,dec from az,el
        tint = rs.tint
        lat = rs.gallat
        total = total + tint
        # now add to histogram 
        glh.add_measurement( lat, tint)
        if 50 * int(iii/50) == iii:
            print("File %5d: az,el = %.3f,%.3f => %.3f,%.3f" % (iii, rs.telaz, rs.telel, rs.gallon, rs.gallat))
        iii = iii + 1

    glh.plot()
    return

if __name__ == "__main__":
    main()



    
