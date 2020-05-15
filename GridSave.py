"""
Model to use the GridClass to make an image of radio astronomical observations
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 20MAY15 GIL one more time on coordiante header updates in FITS
# 20MAY14 GIL update coordinates in FITS header
# 20MAY13 GIL update coordinates in FITS header
# 20MAY12 GIL compute gain ratios for different images.
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
from astropy.wcs import wcs
import radioastronomy
import gainfactor

EPSILON = 0.01
# minimum intensity to include for gain normalization
EPSILON = 1000.0
doRatio = False

# special telescope factors for 2020 March + April 
telescopefactors = [ 1.05216, 0.94350, 1.02153, 0.98935]
telescopefactors = [ 1.09216, 0.94350, 1.02153, 0.98935]
tsum = 0.
for iii in range( 4):
    tsum = tsum + telescopefactors[iii]
tsum = tsum/4.

print( "Check of Telescope Factors: TSum = %f" % (tsum))

def fixImageCoordinates( filename, projection):
    """
    fixImageCoordinates() interpolates pixels to proper reference frame
    """
    printcount = 0

    inname = filename
    nchar = len(inname)
    # strip off the end of the previous image and add the new projection name
    outname = inname[0:nchar-8]
    outname = outname + projection + ".fit"

# get the input image coordinate transform, Usually Cartesian
    win = wcs.WCS(filename)

    hdu = fits.open(filename)[0]
    imageData = fits.getdata( filename)
    imageCopy = copy.deepcopy( imageData)
#
    header = hdu.header
    nx = header['NAXIS1']
    ny = header['NAXIS2']

    crval1 = header['CRVAL1']
    crval2 = header['CRVAL2']
    crpix1 = header['CRPIX1']
    crpix2 = header['CRPIX2']
    cdelt1 = header['CDELT1']
    cdelt2 = header['CDELT2']
    ctype1 = header['CTYPE1']
    ctype2 = header['CTYPE2']

    xmin = crval1 + (1. - crpix1)*cdelt1
    xmax = crval1 + (nx - crpix1)*cdelt1
    ymin = crval2 + (1. - crpix2)*cdelt2
    ymax = crval2 + (nx - crpix2)*cdelt2

    print( "fixImage: %.2f,%2f %.1f,%.1f %.3f,%.3f" % (crval1,crval2,crpix1,crpix2,cdelt1,cdelt2))
    print( "fixImage: %s,%s" % (ctype1,ctype2))
    # redefine for new projection desired
    ctype1 = ctype1[0:4]+projection
    ctype2 = ctype2[0:4]+projection
    print( "=> %s, %s" % (ctype1, ctype2))

    header['CTYPE1'] = ctype1
    header['CTYPE2'] = ctype2

# for output image the reference coordinate x pixel can be anywhere
#    move the center to zero 
    header['CRVAL1'] = 0.

    header['LONPOLE'] = 0.0
    header['LATPOLE'] = 90.0
    header.update()

    tempname = "GridSave.fits"
    hdu = fits.PrimaryHDU(header=header, data=imageCopy)
    print("Outname: %s" % (tempname))
    if os.path.exists(tempname):
        os.remove(tempname)
    hdu.writeto(tempname)

    wout = wcs.WCS(tempname)
    # now that coordinates are defined, remove temporary file
    if os.path.exists(tempname):
        os.remove(tempname)

    pixin = np.array([[0, 0], [nx-1, ny-1]], dtype=np.float64)
    pixout = np.array([[0, 0], [nx-1, ny-1]], dtype=np.float64)

    print("NX, NY: %d,%d" % (nx, ny))

    nan = float("NAN")
#    print("Nan = %f" % (nan))
# assume no data until found
    for jjj in range (ny):
        for iii in range (nx):
            imageCopy[jjj][iii] = nan

# now for output image check all pixel values
    for jjj in range (ny):
        for iii in range (nx):
            # if this image pixal has no value
           pixout[0] = (iii,jjj)
           oworld = wout.wcs_pix2world(pixout, 0)
           xy = oworld[0]
           if np.isnan(xy[0]):
               continue
#                print("pixout: %d,%d : world %.f,%.2f" % (iii,jjj,xy[0],xy[1]))
           pixin[0] = oworld[0]
           ipixels = win.wcs_world2pix(pixin, 0)
# get input pixels for coordinate
           ixy = ipixels[0]
# if outside of current image skip this pixel
           if np.isnan( ixy[0]):
               continue
           ix = int(ixy[0])
           iy = int(ixy[1])
           ix = max( min( nx-1, ix), 0)
           iy = max( min( ny-1, iy), 0)
           ix = int(ix)
           iy = int(iy)
#                print("pixin : %d,%d : world %.f,%.2f" % (ix,iy,xy[0],xy[1]))
#            print("OX,OY:%d,%d <= IX,IY:%d,%d" %( ox,oy, ix,iy))
           imageCopy[jjj][iii] = imageData[iy][ix]

    print("Preparing to write new coordiante transform: %s" % (outname))
    if os.path.exists(outname):
        os.remove(outname)
    newhdu = fits.PrimaryHDU(header=header, data=imageCopy)
    newhdu.writeto(outname)
    print("Wrote new")

    return
    
def writeFitsImage( rs, cpuIndex, grid, projection):
    """
    writeFitsImage() takes a spectrum for describing the observation and a 2 dimensinoal
    array of image data and writes a FITS image
    This program produces two images.   It expects an grid that is in cartisian format.
    The second format described by the input: projection
    """

#    print("Image: ", imageData)
    
    imageData = grid.image
    size = imageData.shape
    imageCopy = copy.deepcopy( imageData)
    nx = size[1]
    ny = size[0]

    # now flip the Y axis of the image to match the FITS Convention
    iy = ny - 1
    for iii in range(ny):
        imageCopy[iii][:] = imageData[iy][:]
        iy = iy - 1

    pixcrd = np.array([[0, 0], [24, 38]], dtype=np.float64)

    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = wcs.WCS(naxis=2)

    gridtype = grid.gridtype.upper()
    print("Grid Type: %s %d" % (gridtype, gridtype.find('RA')))
#    gridtype = "RA"
    if gridtype.find('RA') > -1:
        maptype = 'RA'
        XTYPE = 'RA--'
        YTYPE = 'DEC-'
    else:
        maptype = 'GAL'
        XTYPE = 'GLON'
        YTYPE = 'GLAT'
    xstart = 360.
    ystart = 90.

# select the projection here:
#    projection = "-CYP"
#    projection = "-CAR"

    crval1 = grid.crval1
    crval2 = grid.crval2
    crpix1 = grid.crpix1
    crpix2 = grid.crpix2
    cdelt1 = grid.cdelt1
    cdelt2 = grid.cdelt2
    print('--------- Grid Type: %s (%f,%f %f,%f ' % (gridtype, crval1, crval2, cdelt1, cdelt2))

    hdu = fits.PrimaryHDU()
    header = hdu.header

    dateobs = "%s" % (rs.utc)
    dateobs = dateobs.replace(" ","T")
    mydate = datetime.datetime.now()
    mydate = "%s" % (mydate)
    mydate = mydate[2:10]
    mydate.replace('-','/')

    header['NAXIS1'] = int(nx)
    header['NAXIS2'] = int(ny)
    header['BUNIT'] = 'K-km/s/BEAM'
    maptype = "RA"
    if maptype[0:2] == "RA":
        maptype = "RA"
        header['CTYPE1'] = 'RA---CAR'
    else:
        maptype = "GAL"
        header['CTYPE1'] = 'GLON-CAR'

    # create a cartesian x centered iamge 
    header['CRPIX1'] = nx/2.
    header['CRVAL1'] = 180.
    grid.crval1 = header['CRVAL1']
    header['CDELT1'] = cdelt1
    header['CUNIT1'] = 'deg'
    header['CRVAL2'] = (grid.ymax+grid.ymin)/2.
    grid.crval2 = header['CRVAL2']
    header['CRPIX2'] = ny/2.
    header['CDELT2'] = cdelt2
    header['CUNIT2'] = 'deg'

    grid.gridtype = maptype
    if maptype[0:2] == "RA":
        print("RA: writeFits: %s" % (maptype))
        header['CTYPE2'] = 'DEC--CAR'
    else:
        print("GAL: writeFits: %s" % (maptype))
        header['CTYPE2'] = 'GLAT-CAR'

    header['WCAXES'] = 2
    header['RADESYS'] ='FK5'

# temporarily replace ref coordinate iwth zero
    crval2 = header['CRVAL2']
    crpix2 = header['CRPIX2']
# redefine the reference for the best cartisian format 
    referencevalue = 0.
    dpix = (referencevalue - crval2)/cdelt2
    crpix2 = crpix2 + dpix
# change x axis
    header['CRVAL2'] = referencevalue
    header['CRPIX2'] = crpix2

    header['EQUINOX'] = 2.000000000000E+03 # Equinox of equatorial coordinates
    header['BMAJ'] = 18.1 # Beam major axis in degrees: 80cm horn at 21.1cm
    header['BMIN'] = 18.1 # Beam minor axis in degrees
    header['BPA'] = 0.000000000000E+00 # Beam position angle in degrees
    header['RESTFRQ'] = 1.42040575177E+09 # Line rest frequency, Hz
    header['RESTWAV'] = 0.211061140551    # Line wavelength (m)
    header['DATE-OBS'] = dateobs
    header['DATE'] = mydate
    header['OBSERVER'] = 'Science Aficionado'
    header['OBJECT'] = 'Milky Way'
    header['TELESCOP'] = 'Aficionado Horn'
    header['HISTORY'] = "GridSave.py -- Glen Langston -- 20 May 13"
    header['HISTORY'] = "Observations in March + April 2020"

#    while len(header) < (36 * 4 - 1):
#        header.append()  # Adds a blank card to the end
#    header.delval("EXTEND")
    header.update()

#    hdu = fits.PrimaryHDU(header=header, data=imageData)
    hdu = fits.PrimaryHDU(header=header, data=imageCopy)

    # As file at filePath is deleted now, so we should check if file exists or not not before deleting them
    outname = ("Aficionado_T%d" % (cpuIndex)) + "-" + maptype + projection + ".fit"
    if os.path.exists(outname):
        os.remove(outname)
    hdu.writeto(outname)

# create a second file with new projection
    fixImageCoordinates( outname, projection)

    return

def gridratio( grid1, grid2):
    """
    gridratio computes the ratio of two grids when the values in both grids are non-zero
    This function is used to compute gain ratios
    The average and rms of the ratios are provided along as the grid of ratios
    """

    nx1 = grid1.img_width
    ny1 = grid1.img_height
    nx2 = grid2.img_width
    ny2 = grid2.img_height

    ratio = 0.
    rms = 0.

    if nx1 != nx2:
        print("GridRatio: Nx1 != Nx2 (%d, %d)" % (nx1, nx2))
        return ratio, rms

    if ny1 != ny2:
        print("GridRatio: Ny1 != Ny2 (%d, %d)" % (ny1, ny2))
        return ratio, rms

    count = 0
    nonzero = np.zeros(nx1*ny1)

    # copy to ratio array
    gridratio = copy.deepcopy( grid1)

    for iii in range(nx1):
        for jjj in range(ny1):
            # put in zero as default
            gridratio.image[jjj,iii] = 0.
            if grid1.image[jjj,iii] > EPSILON:
                if grid2.image[jjj,iii] > EPSILON:
                    nonzero[count] = grid1.image[jjj,iii]/grid2.image[jjj,iii]
                    count = count + 1
    if count < 2:
        print ("No overlap in non-zero samples")
        return ratio, rms, gridratio

    nonzero = nonzero[0:count]
    asum = np.sum( nonzero)
    ratio = asum/float(count)
    rms = np.std( nonzero)
    print ("Grid Ratio: %.4f +/- %.4f for %d samples" % (ratio, rms/np.sqrt(count), count))
    # return the ratio grid 
    return ratio, rms, gridratio

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
    FWHM = 3.0  # degrees
    FWHM = 1.0  # degrees
    weight = 1.

    nargs = len(sys.argv)
    if nargs < 2:
        print('GR: GRid Observations of integrated intensity produced by the T Command')
        print('GR produces fits images for each of the horns used for the observations.')
        print('For observations at the same coordinates, the ratios of intensities are also produced.')
        print('The FITS format files require header information, which is copied from the')
        print('Cold Load File provided by the user')
        print('GR RA|GAL <cold file name> <savefile1> [<savefile2> ... <savefileN>]')
        print("")
        print('Glen Langston, National Science Foundation -- 20 May 12')
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
        maptype = 'RA'
    elif gridtype == '-RA':
        xmin = 0.
        xmax = 360.
        ymin = -40.
        ymax = 90.
        xsign = -1.
        xoffset = 360.  # when x = 360. should be at zero.
        maptype = 'RA'
    elif gridtype == '-EL':
        xmin = 0.
        xmax = 360.
        ymin = 0.
        ymax = 90.
        xsign = -1.
        xoffset = 360.  # when x = 360. should be at zero.
        maptype = 'AZEL'
    elif gridtype == 'RA0':
        xmin = 0.
        xmax = 360.
        ymin = -41.
        ymax = 89.
        xsign = -1.
        xoffset = 180.  # when x = 360. should be at zero.
        gridtype = 'RA'
    elif gridtype == 'GAL':
        xmin = -180.
        xmax = 180.
        ymin = -90.
        ymax = 90.
        maptype = 'GAL'

    if gridtype != 'RA' and gridtype != 'GAL' and gridtype != '-RA' and gridtype != "RA0":
        print('Error parsing grid type: ', gridtype)
        print('1st argument should be either RA, -RA or GAL')
        exit()

    rs = radioastronomy.Spectrum()

    if doRatio: 
    #create the grid with map parameters
        grid1 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="-CAR", gridtype=maptype)
        grid2 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="-CAR", gridtype=maptype)
        grid3 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="-CAR", gridtype=maptype)
        grid4 = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="-CAR", gridtype=maptype)
    # put each telescope in a different grid
        grids = [grid1, grid2, grid3, grid4]

    gridall = GridClass.Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width=width, \
                                height=height, dpi=dpi, FWHM=FWHM, \
                                projection="-CAR", gridtype=maptype)
    

    projection = "-AIT"
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
            gainCorr = telescopefactors[iGrid]
            tsum = tsum * gainCorr
            if gridtype == 'RA':
                if doRatio:
                    grids[iGrid].convolve(ra, dec, tsum, weight)
                gridall.convolve( ra, dec, tsum, weight)
            elif gridtype == '-RA':
                x = (ra*xsign) + xoffset
                if doRatio:
                    grids[iGrid].convolve(x, dec, tsum, weight)
                gridall.convolve( x, dec, tsum, weight)
            elif gridtype == 'RA0':
                x = (ra*xsign) + xoffset
                if x < 0:
                    x = x + xmax
                elif x > xmax:
                    x = x - xmax
                if doRatio:
                    grids[iGrid].convolve(x, dec, tsum, weight)
                gridall.convolve( x, dec, tsum, weight)
            else:
                if doRatio:
                    grids[iGrid].convolve(lon, lat, tsum, weight)
                gridall.convolve( lon, lat, tsum, weight)

            if count == 0:
                print('Convolving Coordinates: ', ra, dec, lon, lat)
                print('Convolving Intensities: ', tsum, tsdv, vave, vsdv)
                print('Convolvign Parameters : ', n, time)
            count = count + 1
        # end reading all lines in save file
        f.close()

    # normalize each of the gridded images
    if doRatio:
        grids[0].normalize()
        grids[1].normalize()
        grids[2].normalize()
        grids[3].normalize()
    gridall.normalize()
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

        plt.imshow(gridall.image, interpolation='nearest', cmap=plt.get_cmap('jet'))

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

    crval2 = (xmin + xmax)/2.
    crval1 = (ymin + ymax)/2.
    cdelt1 = (-1./float(dpi)) - .001
    cdelt2 = (1./float(dpi)) + .001
    if doRatio:
# now show eacsh of the images
        for iGrid in range(4):
            imagetemp = copy.deepcopy(grids[iGrid].image)
            imagetemp2 = copy.deepcopy(grids[iGrid].image)
            kkk = myheight - 1
            for jjj in range(myheight):
                imagetemp[:][kkk] = imagetemp2[:][jjj]
                kkk = kkk - 1
                grids[iGrid].image = imagetemp
            writeFitsImage( rs, iGrid+2, grids[iGrid], projection)

        # put each telescope in a different grid
        ratio1 = copy.deepcopy(grid1)
        ratio2 = copy.deepcopy(grid1)
        ratio3 = copy.deepcopy(grid1)
        gratios = [ratio1, ratio2, ratio3]
        ratios = np.zeros(3)
        rmss = np.zeros(3)

        jGrid = 3
        for iGrid in range(3):
            print("Gain Ratios for Telescopes T%d and T%d" % (iGrid+2, jGrid+2))
            ratio, rms, aratio = gridratio(grids[iGrid], grids[jGrid])
            ratios[iGrid] = ratio
            rmss[iGrid] = rms
            writeFitsImage( rs, iGrid+2, aratio, projection)
    
    writeFitsImage( rs, 0, gridall, projection)
    plt.show()

if __name__ == "__main__":
    main()

#SIMPLE  =                    T / conforms to FITS standard                      
#BITPIX  =                  -32 / array data type                                
#NAXIS   =                    2 / number of array dimensions                     
#NAXIS1  =                 4323                                                  
#NAXIS2  =                 2163                                                  
#OBJECT  = 'HI4PI   '           / The HI 4-PI Survey                             
#TELESCOP= 'Effelsberg 100m RT; ATNF Parkes 64-m' / Telescope names              
#ORIGIN  = 'AIfA/MPIfR Bonn; ATNF Sydney' / Organisations or Institutions        
#REFERENC= 'HI4PI Collaboration 2016' / A&A                                      
#RESTFRQ =        1420405751.77                                                  
#RESTWAV =       0.211061140541                                                  
#CDELT1  = -0.08333333330000001                                                  
#CRPIX1  =               2162.0                                                  
#CRVAL1  =                  0.0                                                  
#CTYPE1  = 'RA---CAR'                                                            
#CUNIT1  = 'deg     '                                                            
#CDELT2  =  0.08333333330000001                                                  
#CRPIX2  =               1082.0                                                  
#CRVAL2  =                  0.0                                                  
#CTYPE2  = 'DEC--CAR'                                                            
#CUNIT2  = 'deg     '                                                            
#WCSAXES =                    2                                                  
#RADESYS = 'FK5     '                                                            
#EQUINOX =               2000.0                                                  
#LONPOLE =                  0.0                                                  
#LATPOLE =                 90.0                                                  
#BUNIT   = 'cm^(-2) '                                                            
#BPA     =                  0.0                                                  
#BMAJ    =               0.2706                                                  
#BMIN    =               0.2706                                                  
#VMIN    = -1.5169915133210E+21                                                  
#VMAX    = 2.39529868415209E+22                                                  
#CHECKSUM= '9eqZJcnY9cnYGcnY'   / HDU checksum updated 2016-09-15T23:38:45       
#DATASUM = '3638685465'         / data unit checksum updated 2016-09-15T23:38:45 
#END                                                                             
