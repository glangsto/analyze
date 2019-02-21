"""
Model to use the GridClass to make an image of radio astronomical observations
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 18FEB27 GIL initial version to show position velocity plot
# 18FEB19 GIL possibly grid max
# 17FEB03 GIL add comments and cleanup
# 17JAN28 GIL finish initial version
# 17JAN09 GIL initial version based on mandelbrot() python example

import sys
import numpy as np
#from numba import jit
from matplotlib import pyplot as plt
#from matplotlib import colors
import radioastronomy
import GridClass

def main():
    """
    Main executable for gridding astronomical data
    """
    dpi = 1
    dpi = 2
    width = int(360)
    height = int(180)
    mywidth = int(width*dpi)
    myheight = int(height*dpi)
    FWHM = 1.0  # degrees
    weight = 1.

    nargs = len(sys.argv)
    if nargs < 2:
        print 'GO: Grid Observations'
        print 'GO RA|GAL|VLAT|VTIME datafiles'
        exit()

    gridtype = sys.argv[1]
    gridtype = gridtype.upper()
    print 'Grid Type: ', gridtype

    if gridtype == 'VTIME':
        xmin = 0.
        xmax = 360..
        ymin = -120.
        ymax = 120.

    if gridtype == 'VLON':
        xmin = 0.
        xmax = 360..
        ymin = -120.
        ymax = 120.

    if gridtype == 'VLAT':
        xmin = 0.
        xmax = 360..
        ymin = -120.
        ymax = 120.

    if gridtype != 'VLON' and gridtype != 'VLAT' and gridtype != 'VTIME':
        print 'Error parsing grid type: ', gridtype
        print '1st argument should be either VLON, VLAT or VTIME'
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
    rs = radioastronomy.Spectrum()

    for filename in names:
        
        rs.read_spec_ast(filename)
        rs.azel2radec()

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
            else:
                mygrid.convolve(lon, lat, tsum, weight)
            if count == 0:
                print 'Convolving Coordinates: ', ra, dec, lon, lat
                print 'Convolving Intensities: ', tsum, tsdv, vave, vsdv
                print 'Convolvign Parameters : ', n, time
            count = count + 1

    mygrid.normalize()
#    mygrid.check()
    zmin = 00.
    zmax = 2000.
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
            plt.title("Citizen Science: Observing our Galaxy - %s to %s" % (firsttime, lasttime))
        else:
            plt.title("Citizen Science: Observing our Galaxy - %s" % (firsttime))
        if gridtype == 'RA':
            plt.xlabel("Right Ascension (hours)")
            plt.ylabel("Declination (degrees)")
            labels = ticks/(mywidth/24)
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

