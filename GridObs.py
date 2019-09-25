"""
Model to use the GridClass to make an image of radio astronomical observations
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 19SEP23 GIL slightly bigger map in declination
# 18FEB19 GIL possibly grid max
# 17FEB03 GIL add comments and cleanup
# 17JAN28 GIL finish initial version
# 17JAN09 GIL initial version based on mandelbrot() python example

import sys
import numpy as np
#from numba import jit
from matplotlib import pyplot as plt
#from matplotlib import colors
import GridClass

def main():
    """
    Main executable for gridding astronomical data
    """
    dpi = 1
    dpi = 2
    width = int(360)
    height = int(130)
    mywidth = int(width*dpi)
    myheight = int(height*dpi)
    FWHM = 7.5  # degrees
    FWHM = 1.0  # degrees
    weight = 1.

    nargs = len(sys.argv)
    if nargs < 2:
        print 'GO: Grid Observations'
        print 'GO RA|GAL datafiles'
        exit()

    gridtype = sys.argv[1]
    gridtype = gridtype.upper()
    print 'Grid Type: ', gridtype

    # enable having ra going from 24 to 0 hours == 360 to 0 degrees
    xsign = 1.
    xoffset = 0.
    if gridtype == 'RA':
        xmin = 0.
        xmax = 360.
        ymin = -40.
        ymax = 90.
    elif gridtype == '-RA':
        xmin = 0.
        xmax = 360.
        ymin = -40.
        ymax = 90.
        xsign = -1.
        xoffset = 360.  # when x = 360. should be at zero.
    elif gridtype == 'RA0':
        xmin = 0.
        xmax = 360.
        ymin = -41.
        ymax = 89.
        xsign = -1.
        xoffset = 180.  # when x = 360. should be at zero.
    elif gridtype == 'GAL':
        xmin = -180.
        xmax = 180.
        ymin = -90.
        ymax = 90.

    if gridtype != 'RA' and gridtype != 'GAL' and gridtype != '-RA' and gridtype != "RA0":
        print 'Error parsing grid type: ', gridtype
        print '1st argument should be either RA, -RA or GAL'
        exit()


    #create the grid with map parameters
    mygrid = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="Mercator", gridtype=gridtype)

# first read through all data and find hot load
    names = sys.argv[2:]
    names = sorted(names)

    firsttime = ""
    lasttime = ""
    count = 0
    for filename in names:
        f = open(filename)
        lines = f.readlines()
        f.close()

        for line in lines:
            line = line.strip()
            count = count + 1
            if len(line) < 10:
                continue
            if line[0] == '#':
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                float(parts[0])
            except ValueError:
                print "Not a float"
                print 'Parts: ', parts
                continue
            ra = float(parts[0])
            dec = float(parts[1])
            lon = float(parts[2])
            lat = float(parts[3])
            tsum = float(parts[4])
            tsdv = float(parts[5])
            tmax = float(parts[6])
            vave = float(parts[7])
            vsdv = float(parts[8])
            el = float(parts[9])
            n = float(parts[10])
            time = parts[11]
            if firsttime == "":
                firsttime = time
            else:
                lasttime = time

#            if vave > -100. and vave < 100:
#                mygrid.convolve( lon, lat, vave, 1.)
            if gridtype == 'RA':
                mygrid.convolve(ra, dec, tsum, weight)
            elif gridtype == '-RA':
                x = (ra*xsign) + xoffset
                mygrid.convolve(x, dec, tsum, weight)
            elif gridtype == 'RA0':
                x = (ra*xsign) + xoffset
                if x < 0:
                    x = x + xmax
                elif x > xmax:
                    x = x - xmax
                mygrid.convolve(x, dec, tsum, weight)
            else:
                mygrid.convolve(lon, lat, tsum, weight)
            if count == 0:
                print 'Convolving Coordinates: ', ra, dec, lon, lat
                print 'Convolving Intensities: ', tsum, tsdv, vave, vsdv
                print 'Convolvign Parameters : ', n, time
            count = count + 1

    mygrid.normalize()
#    mygrid.check()
    zmin = -1500.
    zmax = 6000.
# limit grid intensities for plotting
    mygrid.limit(zmin, zmax)

    subplots = False

    if subplots:
        fig, ax = plt.subplots(figsize=(myheight, mywidth), dpi=dpi)

        if gridtype == 'RA':
            cax = fig.add_axes([-180, 180], [-90, 90])
        else:
            cax = fig.add_axes([0, 24], [-90, 90])

#        im = ax.imshow(mygrid.image, interpolation='nearest')
        cbar = fig.colorbar(cax, ticks=[zmin, zmax], orientation='horizontal')
        cbar.ax.set_yticklabels([str(zmin), str(zmax)])

        ax.set_title("Citizen Science: Horn observations of our Galaxy")
    else:
#y_ticks = ymin + (ymax-ymin)*ticks/myheight

        ticks = np.arange(0, mywidth, 30*dpi)
        x_ticks = xmin + ((xmax-xmin)*ticks/mywidth)

        plt.imshow(mygrid.image, interpolation='nearest', cmap=plt.get_cmap('jet'))

        if firsttime != lasttime:
            plt.title("Citizen Science: Observing our Galaxy: %s to %s" % (firsttime, lasttime))
        else:
            plt.title("Citizen Science: Observing our Galaxy: %s" % (firsttime))
        if gridtype == 'RA':
            plt.xlabel("Right Ascension (hours)")
            plt.ylabel("Declination (degrees)")
            labels = ticks/(mywidth/24)
            yticks = np.arange(0, myheight, 15*dpi)
        elif gridtype == '-RA':
            plt.xlabel("Right Ascension (hours)")
            plt.ylabel("Declination (degrees)")
            labels = 24 - (ticks/(mywidth/24))
            labels[0] = 0
            labels[0] = 24
            yticks = np.arange(0, myheight, 15*dpi)
        elif gridtype == 'RA0': # put 0 hours in middle of plot
            plt.xlabel("Right Ascension (hours)")
            plt.ylabel("Declination (degrees)")
            labels = 12 - (ticks/(mywidth/24))
            nlabels = len(labels)
            for iii in range(nlabels):
                if labels[iii] < 0:
                    labels[iii] = 24 + labels[iii]
                if labels[iii] == 24:
                    labels[iii] = 0
            yticks = np.arange(0, myheight, 15*dpi)
        else:
            yticks = np.arange(0, myheight, 30*dpi)
            ticks = np.arange(0, mywidth, 30*dpi)
            x_ticks = xmin + (xmax-xmin)*ticks/mywidth
            labels = x_ticks
            plt.xlabel("Galactic Longitude (degrees)")
            plt.ylabel("Galactic Latitude (degrees)")
        # wnat an integer list of labels
#        slabels = str(labels)
        print ticks, labels
        y_ticks = ymax - (ymax-ymin)*yticks/myheight
        plt.yticks(yticks, y_ticks)
        plt.xticks(ticks, labels, rotation='horizontal')
        plt.colorbar()

    plt.show()

if __name__ == "__main__":
    main()

