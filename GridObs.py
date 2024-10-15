"""
Model to use the GridClass to make an image of radio astronomical observations
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 21JAN06 GIL change use of split
# 19SEP23 GIL slightly bigger map in declination
# 18FEB19 GIL possibly grid max
# 17FEB03 GIL add comments and cleanup
# 17JAN28 GIL finish initial version
# 17JAN09 GIL initial version based on mandelbrot() python example

import sys
import numpy as np
#from numba import jit
from matplotlib import pyplot as plt
import datetime

#from matplotlib import colors
import GridClass
import radioastronomy

# below are the indices into a data line
IDATE = 0
ITIME = 1
ITEL = 2
IAZ = 3
IEL = 4
ITSYS = 5
ITRX = 6
ITRMS = 7
IAVET = 8
IKCOUNT = 9
IPEAKK = 10
IPEAKV = 11
IPEAKVRMS = 12
IVSUM = 13
IVSUMRMS = 14
IINTK = 15
IINTKRMS = 16
ISCALE = 17
IRA = 18
IDEC = 19
IGLON = 20
IGLAT = 21

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
    FWHM = 1.0  # degrees
    FWHM = 10.0  # degrees
    FWHM = 1.0  # degrees
    weight = 1.
    lastaz = -1
    lastel = -1
    
    nargs = len(sys.argv)
    if nargs < 2:
        print('GO: Grid Observations')
        print('GO RA|GAL|EL Spectrum Summary-files')
        exit()

    gridtype = sys.argv[1]
    gridtype = gridtype.upper()
    print('Grid Type: %s' % (gridtype))
    spectrum = sys.argv[2]
    specparts = spectrum.split(".")
    nparts = len(specparts)
    if nparts < 2:
        print("Invalid Spectrum File name: %s" % (spectrum))
        exit()
    else:
        filetype = specparts[nparts-1]
        if filetype != "ast" and filetype != "hot" and filetype != "kel":
            print("Spectrum file type not recognized: %s (%s)" % \
                  (spectrum, filetype))
            exit()   

    # create a radio astronomy spectrum structure to use utilities
    rs = radioastronomy.Spectrum()
    # read the file describing the observation and telescope location
    rs.read_spec_ast(spectrum)
    
    if rs.nChan <= 0:
        print("File %s is not an astronomy, hot or note file" % (spectrum))
        print("The first file sets the telescope location for mapping")
        exit()
        
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
    elif gridtype == '-EL'  or gridtype == 'EL':
        xmin = 0.
        xmax = 360.
        ymin = 0.
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

    if gridtype != 'RA' and gridtype != 'GAL' and gridtype != '-RA' and gridtype != "RA0" and gridtype != 'EL' and gridtype != '-EL':
        print('Error parsing grid type: %s' % (gridtype))
        print('1st argument should be either RA, -RA, RA0, EL, -EL or GAL')
        exit()

    #create the grid with map parameters
    mygrid = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="Mercator", gridtype=gridtype)
    # use integer to match gridding type
    gridindex = mygrid.gridindex
# setup to parse the date and time
    timefmt = "%Y-%m-%d %H:%M:%S"

    print("Grid Class: %s == %d " % (gridtype, gridindex))
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

        count = 0
        for line in lines:
            line = line.strip()
            if len(line) < 10:
                continue
            if line[0] == '#':
                continue
            parts = line.split()
            if len(parts) < 5:
                continue

            # test that we partially understand the data format
            try:
                float(parts[IAZ])
            except ValueError:
                print("Azimuth Not a float")
                print(('Parts: ', parts[0:IAZ+1]))
                continue
            # Transfer
            date = "20"+parts[IDATE]
            time = parts[ITIME]
#            timefmt = "%Y-%m-%d %H:%M:%S"
            utc = datetime.datetime.strptime(date + " " + time, timefmt)
#            utc = datetime.fromisoformat(date + " " + time)
            rs.utc = utc

            telIndex = int(parts[ITEL])
            rs.telaz = float(parts[IAZ])
            rs.telel = float(parts[IEL])
            if count == 0:
                lastaz = rs.telaz
                lastel = rs.telel
                print("First Az,el = %.2f, %.2f " % (rs.telaz, rs.telel))
            else:
                if rs.telaz != lastaz or rs.telel != lastel:
                    print("New Az,el = %.2f, %.2f " % (rs.telaz, rs.telel))
                    lastaz = rs.telaz
                    lastel = rs.telel
                    
            # now recompute ra,dec
            rs.azel2radec()
            ra = rs.ra
            dec = rs.dec

            tsys = float(parts[ITSYS])
            trx = float(parts[ITRX])
            trms = float( parts[ITRMS])
            avet = int(parts[IAVET])

            lon = rs.gallon
            lat = rs.gallat
            tsum = float(parts[IINTK])
            tsdv = float(parts[IINTKRMS])
            tmax = float(parts[IPEAKK])
            vave = float(parts[IPEAKV])
            vsdv = float(parts[IPEAKVRMS])
            n = float(parts[IAVET])
            time = parts[0]
#            time = parts[11]
            if firsttime == "":
                firsttime = time
            else:
                lasttime = time

            if gridindex == GridClass.GRIDRADEC:
                mygrid.convolve(ra, dec, tsum, weight)
            elif gridindex == GridClass.GRIDMRADEC:
                x = (ra*xsign) + xoffset
                mygrid.convolve(x, dec, tsum, weight)
            elif gridindex == GridClass.GRIDEL:
                x = ra + xoffset
                if x < 0:
                    x = x + xmax
                elif x > xmax:
                    x = x - xmax
                mygrid.convolve(x, rs.telel, tsum, weight)
            elif gridtype == GridClass.GRID0RADEC:
                x = (ra*xsign) + xoffset
                if x < 0:
                    x = x + xmax
                elif x > xmax:
                    x = x - xmax
                mygrid.convolve(x, dec, tsum, weight)
            else:  # gridding galactic cordinates
                mygrid.convolve(lon, lat, tsum, weight)
            # Show progress
            if count%300 == 0:
                print('Convolving Coordinates: %6.2f,%6.2f => %6.2f,%6.2f' % \
                       (ra, dec, lon, lat))
                print('Convolving Intensities: %6.2f,%6.2f' % (tsum, tsdv))
                print('Convolvign Parameters : %s' % (time))

            count = count + 1

    if count < 1:
        print("No samples gridded, Exiting")
        exit()

    mygrid.normalize()
    mygrid.check()
    zmin = -500.
    zmax = 1000.
# limit grid intensities for plotting
    mygrid.set_ij( 0, 0, zmax, 1.)
    mygrid.set_ij( 1, 1, zmin, 1.)
#    mygrid.limit(zmin, zmax)

    subplots = False

    if subplots:
        fig, ax = plt.subplots(figsize=(myheight, mywidth), dpi=dpi)

        if gridtype == 'RA':
            cax = fig.add_axes([-180, 180], [-90, 90])
        elif gridtype == 'EL'  or gridtype == '-EL':
            cax = fig.add_axes([0, 24], [0, 90])
        else:
            cax = fig.add_axes([0, 24], [-90, 90])

#        im = ax.imshow(mygrid.image, interpolation='nearest')
        cbar = fig.colorbar(cax, ticks=[zmin, zmax], \
                            orientation='horizontal', label="K/km/s")
        cbar.ax.set_yticklabels([str(zmin), str(zmax)])

        ax.set_title("Citizen Science: Horn observations of our Galaxy")
    else:
#        y_ticks = ymin + (ymax-ymin)*ticks/myheight

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
        elif gridtype == '-EL' or gridType == 'EL':
            plt.xlabel("Right Ascension (hours)")
            plt.ylabel("Elevation (degrees)")
            labels = 24 - (ticks/(mywidth/24))
            labels[0] = 0
            labels[0] = 24
            yticks = np.arange(0, 90, 10)
#            yticks = np.arange(0, int(myheight/9), 30*dpi)
#            yticks = np.arange(0, myheight, 10*dpi)
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
        print((ticks, labels))
        y_ticks = ymax - (ymax-ymin)*yticks/myheight
        if gridtype == 'EL' or gridtype == '-EL':
            print((yticks, y_ticks))
            plt.yticks(yticks, y_ticks)
        else:
            plt.yticks(yticks, y_ticks)
        plt.xticks(ticks, labels, rotation='horizontal')
        plt.colorbar()

    plt.show()

if __name__ == "__main__":
    main()
