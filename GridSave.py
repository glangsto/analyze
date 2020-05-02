a"""
Model to use the GridClass to make an image of radio astronomical observations
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 20MAY01 GIL grid different telescopes to different images for gain and offset merging.
# 20APR30 GIL initial version based on GridObs.py

import sys
import os
import copy
import numpy as np
from matplotlib import pyplot as plt
import datetime
import GridClass
from astropy.io import fits
import radioastronomy
import gainfactor

def writeFitsImage( rs, cpuIndex, crval1, crval2, cdelt1, cdelt2, imageData):
    """
    writeFitsImage() takes a spectrum for describing the observation and a 2 dimensinoal
    array of image data and writes a FITS image
    """

#    print("Image: ", imageData)
    
    size = imageData.shape
    print("Data size: ", size)
    nx = size[0]
    ny = size[1]

    hdu = fits.PrimaryHDU(imageData)

    header = hdu.header
    while len(header) < (36 * 4 - 1):
        header.append()  # Adds a blank card to the end

    header['NAXIS1'] = int(size[0])
    header['NAXIS2'] = int(size[1])
    header['CRPIX1'] = size[0]/2.
    header['CRPIX2'] = size[1]/2.
    header['CRVAL1'] = crval1
    header['CRVAL2'] = crval2
    header['CDELT1'] = cdelt1
    header['CDELT2'] = cdelt2
    header['OBSERVER'] = 'Science Aficionado'
    header['OBJECT'] = 'Milky Way'
    
    outname = ("AficionadoMap_T%d" % (cpuIndex)) + ".fit"
    hdu.writeto(outname)
    return

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
    FWHM = 10.0  # degrees
    FWHM = 5.0  # degrees
    FWHM = 1.0  # degrees
    weight = 1.

    nargs = len(sys.argv)
    if nargs < 2:
        print('GO: Grid Observations')
        print('GO RA|GAL datafiles')
        exit()

    gridtype = sys.argv[1]
    gridtype = gridtype.upper()
    print('Grid Type: ', gridtype)

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
    elif gridtype == '-EL':
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

    if gridtype != 'RA' and gridtype != 'GAL' and gridtype != '-RA' and gridtype != "RA0":
        print('Error parsing grid type: ', gridtype)
        print('1st argument should be either RA, -RA or GAL')
        exit()

    rs = radioastronomy.Spectrum()

    #create the grid with map parameters
    grid1 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="Mercator", gridtype=gridtype)
    grid2 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="Mercator", gridtype=gridtype)
    grid3 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="Mercator", gridtype=gridtype)
    grid4 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="Mercator", gridtype=gridtype)
    
    # put each telescope in a different grid
    grids = [grid1, grid2, grid3, grid4]

# coldfile 
    coldfile = sys.argv[2]
# get telescope geographic location etc
    print("Reading Observing parameters from: %s" % (coldfile))
    rs.read_spec_ast(coldfile)
    print("Observer: %s " % (rs.observer))

# first read through all data and find hot load
    names = sys.argv[3:]
    names = sorted(names)

    firsttime = ""
    lasttime = ""
    count = 0
    # setup grid indicies so that cpuIndex goes to the correct grid
    # This assumes telescopes 2,3,4,5 are being used] 
    gridIndex = [0,0,0,1,2,3]
    # for all save Files to Grid
    for filename in names:
        print("File: %s" % (filename))
        f = open(filename)

        date = "Unknown"
        while date != "":
            date, time, cpuIndex, telaz, telel, tSys, tRx, tRms, tint, KperC, tSourcemax, velSource, dV, tVSum, tVSumRms, tSumKmSec, dTSumKmSec, gainFactor = gainfactor.readSaveValues( f)
            dlen = len(date)
            if dlen < 1:
                break
            if date[0] == "#":
                continue
            # else not a comment process the line
            count = count + 1
            isodate = "20"+date+"T"+time
#            print("DateTime: %s" % (isodate))
            rs.utc = datetime.datetime.strptime(isodate,"%Y-%m-%dT%H:%M:%S")
#            print("Utc: %s" % (rs.utc))
            rs.telaz = telaz
            rs.telel = telel
            rs.azel2radec()

            ra = rs.ra
            dec = rs.dec
            lon = rs.gallon
            lat = rs.gallat
            tsum = tSumKmSec
            tsdv = dTSumKmSec
            tmax = tSourcemax
            vave = tVSum
            vsdv = tVSumRms
            if firsttime == "":
                firsttime = date
            else:
                lasttime = date

#            if vave > -100. and vave < 100:
#                mygrid.convolve( lon, lat, vave, 1.)
            iGrid = gridIndex[cpuIndex]
            if gridtype == 'RA':
                grids[iGrid].convolve(ra, dec, tsum, weight)
            elif gridtype == '-RA':
                x = (ra*xsign) + xoffset
                grids[iGrid].convolve(x, dec, tsum, weight)
            elif gridtype == 'RA0':
                x = (ra*xsign) + xoffset
                if x < 0:
                    x = x + xmax
                elif x > xmax:
                    x = x - xmax
                grids[iGrid].convolve(x, dec, tsum, weight)
            else:
                grids[iGrid].convolve(lon, lat, tsum, weight)

            if count == 0:
                print('Convolving Coordinates: ', ra, dec, lon, lat)
                print('Convolving Intensities: ', tsum, tsdv, vave, vsdv)
                print('Convolvign Parameters : ', n, time)
            count = count + 1
        # end reading all lines in save file
        f.close()

    # normalize each of the gridded images
    grids[0].normalize()
    grids[1].normalize()
    grids[2].normalize()
    grids[3].normalize()
#    mygrid.check()
#    zmin = -1000.
#    zmax = 3000.
# limit grid intensities for plotting
#    mygrid.set_ij( 0, 0, zmax, 1.)
#    mygrid.set_ij( 1, 1, zmin, 1.)
#    mygrid.limit(zmin, zmax)

    subplots = False

    if subplots:
        fig, ax = plt.subplots(figsize=(myheight, mywidth), dpi=dpi)

        if gridtype == 'RA':
            cax = fig.add_axes([-180, 180], [-90, 90])
        else:
            cax = fig.add_axes([0, 24], [-90, 90])

        cbar = fig.colorbar(cax, ticks=[zmin, zmax], orientation='horizontal')
        cbar.ax.set_yticklabels([str(zmin), str(zmax)])

        ax.set_title("Citizen Science: Horn observations of our Galaxy")
    else:
#y_ticks = ymin + (ymax-ymin)*ticks/myheight

        ticks = np.arange(0, mywidth, 30*dpi)
        x_ticks = xmin + ((xmax-xmin)*ticks/mywidth)

        plt.imshow(grids[0].image, interpolation='nearest', cmap=plt.get_cmap('jet'))

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
        elif gridtype == '-EL':
            plt.xlabel("Right Ascension (hours)")
            plt.ylabel("Elevation (degrees)")
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
        print(ticks, labels)
        y_ticks = ymax - (ymax-ymin)*yticks/myheight
        plt.yticks(yticks, y_ticks)
        plt.xticks(ticks, labels, rotation='horizontal')
        plt.colorbar()


    crval1 = (xmin + xmax)/2.
    crval2 = (ymin + ymax)/2.
    cdelt1 = 1./float(dpi)
    cdelt2 = 1./float(dpi)
    for iGrid in range(4):
        imagetemp = copy.deepcopy(grids[iGrid].image)
        imagetemp2 = copy.deepcopy(grids[iGrid].image)
        kkk = myheight - 1
        for jjj in range(myheight):
            imagetemp[:][kkk] = imagetemp2[:][jjj]
            kkk = kkk - 1
        grids[iGrid].image = imagetemp
        writeFitsImage( rs, iGrid+2, crval1, crval2, cdelt1, cdelt2, grids[iGrid].image)
    plt.show()

if __name__ == "__main__":
    main()

