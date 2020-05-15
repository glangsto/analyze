"""
Grid astronomical observations producing a plot of either RA-DEC
or Galactic Longitude-Latitude
Glen Langston
"""
# Functions to create a grid and place astronomical data on that
# grid with a convolving function
# HISTORY
# 20MAY13 GIL Add DATE-OBS
# 19OCT03 GIL make grid pattern circular
# 17FEB03 GIL add comments and cleanup
# 17JAN28 GIL finish initial version
# 17JAN09 GIL initial version based on mandelbrot() python example

import numpy as np
import angles
from matplotlib import pyplot as plt

#matplotlib inline

class Grid(object):
    """
    Grid image class for placing weighted spectral intensities on a grid.
    This class defines an image associated with the data and a weights image.
    The steps are:
    1. Init the empty images based on input parameters
    2. Compute the transform from angular coordiates to grid location
    3. Compute the weighting function for the coordinate value and grid offset
    4. Convolve the intensity with the weighting function
    5. Sum to grid and weights images
    6. repeat for all input values.
    7. Normalise grid
    8. Compute the time variation of gains (optional)
    """
    def __init__(self, xmin=0., xmax=360., ymin=-90., ymax=90., \
                 width=360, height=180, dpi=1, FWHM=13.0, \
                     projection="-CAR", gridtype="RA"):
        """
        initialize all Grid class values
        many will be overwritten later
        """
        self.img_width = int(width*dpi)
        self.img_height = int(height*dpi)
        self.width = int(width)
        self.height = int(height)
        self.dpi = dpi # pixels per unit: eg 10 = 10 pixels per degree == 0.1 degree per pixel
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        # FITS coordinate info
        self.crval2 = (xmin + xmax)/2.
        self.crval1 = (ymin + ymax)/2.
        self.cdelt1 = (-1./float(dpi)) - .001
        self.cdelt2 = (1./float(dpi)) + .001
        self.crpix1 = width/2.
        self.crpix2 = height/2.
        
        self.xminrad = xmin*np.pi/180.  # convert to radians
        self.xmaxrad = xmax*np.pi/180.
        self.yminrad = ymin*np.pi/180.
        self.ymaxrad = ymax*np.pi/180.
        # prepare to compute pixel coordinate from x, y inputs
        self.dx = (self.xmax - self.xmin)
        self.dy = (self.ymax - self.ymin)
        self.idx = float(self.img_width) / self.dx
        self.idy = float(self.img_height) / self.dy
        self.image = np.empty((self.img_height, self.img_width), dtype=np.float32)
        self.weights = np.empty((self.img_height, self.img_width))
        self.FWHM = FWHM # degrees
#        self.sigma = self.FWHM/(2*np.sqrt(2*np.log(2.)))
        self.sigma = self.FWHM/(2.*np.sqrt(2.*np.log(2.)))
        self.sigma2 = self.sigma*self.sigma
#        self.WF = -0.5/self.sigma2  # used in exponent of gaussian
        self.WF = -.1/self.sigma2  # used in exponent of gaussian
        self.gridtype = str(gridtype)  # Either RA or GAL

#        print 'Init: x min,max  (d) : ', self.xmin, self.xmax
#        print 'Init: x min,max (rad): ', self.xminrad, self.xmaxrad
        self.projection = projection
        self.printcount = 0
        # check that the image is initiazed correctly
        self.zero()

    def weightfunction(self, r):
        """
        Returns the gaussian weighting funciton for a complex coordinate distance
        r = Input radius (degrees, or the same units as sqrt(WF)).
        """
        return np.exp(r*r*self.WF)

    def save_image(self, fig):
        """
        Save the gridded image as picture. Work still in progress!
        """
        filename = "grid_%dx%d.png" % (self.width, self.height)
        print(fig)
        f = open(filename, 'w')
        f.close()
#        fig.savefig(filename)

    def ii(self, x):
        """
        return the x index grid coordinate for an input point
        """
        ix = self.idx * (x - self.xmin)  # replace with full coordinate projection
        ix = int(ix)
        if ix < 0:
            ix = ix + self.img_width
        elif ix >= self.img_width:
            ix = ix - self.img_width

#        if x > 350.:
#            print "x,ix", x, ix,'dx: %8.3e, x min,max: %6.1f,%6.1f' % \
#               ( self.idx, self.xmin, self.xmax)

        return ix

    def jj(self, y):
        """
        return the y index grid coordinate for an input point
        """
        iy = int(self.idy * (self.ymax - y))  # preplace with full coordinate calc
        return iy

    def xx(self, ix):
        """
        return the x icoordinate for an input point grid index
        """
        return ((ix/self.idx) + self.xmin)

    def yy(self, iy):
        """
        return the x coordinate for an input point grid index
        """
        return (self.ymax - (iy/self.idy))

    def xy(self, ix, iy, z, r, inweight):
        """
        xy computes the weight for a value x y z at a distance r from the measurement
        x - x-coordinate (degrees)
        y - y-coordinate (degrees)
        z - intensity of a position x',y', a distance r from x,y
        """

        convolutionweight = float(self.weightfunction(r))
        totalweight = float(inweight * convolutionweight)
        if self.printcount == 0:
            print('  i    j    T (K km/s)   R (d)    Weight ')
        if self.printcount < 1 or (not np.isfinite(totalweight)) \
                or (50000*int(self.printcount/50000) == self.printcount) \
                or ((50000*int(self.printcount/50000)+1) == self.printcount):
            print('%4d,%4d: %8.1f %9.2f   %8.3e' % (ix, iy, z, r, totalweight))

        self.printcount = self.printcount + 1

        if ix >= self.img_width:
            ix = ix - self.img_width
        elif ix < 0:
            ix = ix + self.img_width
        if iy >= self.img_height:
            iy = iy - self.img_height
        elif iy < 0:
            iy = iy + self.img_height
        self.image[iy, ix] = self.image[iy, ix] + (totalweight*z)
        self.weights[iy, ix] = self.weights[iy, ix] + totalweight

    def convolve(self, x, y, z, inweight):
        """
        convolve a measurement onto a range of grid coordinates
        """
        nwidth = int(4.1 * self.FWHM * self.dpi)

        ix = self.ii(x)
        jy = self.jj(y)
        x0 = angles.d2r(x) # convert to radians
        if x0 < self.xminrad:
            x0 = x0 + (2.*np.pi)
        elif x0 > self.xmaxrad:
            x0 = x0 - (2.*np.pi)
        y0 = angles.d2r(y)

        for iii in range(-nwidth, nwidth):
            iix = ix + iii
            if iix < 0:
                iix = self.img_width + iix
            elif iix >= self.img_width:
                iix = iix - self.img_width
            xx = self.xx(iix)
            rx = angles.d2r(xx)  # convert to radians

            for jjj in range(-nwidth, nwidth):
                jjy = jy + jjj
                if jjy < 0:
                    jjy = self.img_height + jjy
                elif jjy >= self.img_height:
                    jjy = jjy - self.img_height
                yy = self.yy(jjy)
                ry = angles.d2r(yy) # conver to radians

                if ry < self.yminrad:
                    ry = ry + (2.*np.pi)
                elif ry > self.ymaxrad:
                    ry = ry - (2.*np.pi)

                # finally get angular separation in degrees
                r = angles.r2d(angles.sep(x0, y0, rx, ry))
#                if r > 4.*self.FWHM:  # round convolving function
                if r > 3.*self.FWHM:  # round convolving function
                    continue

                # add the colvolved measurement to the grid.
                self.xy(iix, jjy, z, r, inweight)
        return
                
    def normalize(self):
        """
        normalize() computes the image from the sum and weights images
        """
        nonzero = 0
        for iii in range(self.img_width):
            for jjj in range(self.img_height):
                if self.weights[jjj, iii] != 0.:
                    self.image[jjj, iii] = \
                        self.image[jjj, iii]/self.weights[jjj, iii]
                    nonzero = nonzero + 1
        print('Normalize: ', nonzero, ' Values in the grid')

    def check(self):
        """
        check() checks the values for an image to make sure they are all finite
        """

        nonzero = 0
        for iii in range(self.img_width):
            for jjj in range(self.img_height):
                if not np.isfinite(self.image[jjj, iii]):
                    print('Non-finite value at: ', iii, jjj, self.image[jjj, iii])
                    self.image[jjj, iii] = 0
                    nonzero = nonzero + 1
                if not np.isfinite(self.weights[jjj, iii]):
                    print('Non-finite weight at: ', iii, jjj, self.weights[jjj, iii])
                    self.weights[jjj, iii] = 0
                    nonzero = nonzero + 1
        print('Check: ', nonzero, ' Illegal Values in the grid')

        return

    def limit(self, ymin, ymax):
        """
        check() checks the values for an image to make sure they are all finite
        """

        datamin = self.image[0, 0]
        datamax = self.image[0, 0]

        for iii in range(self.img_width):
            for jjj in range(self.img_height):
                if self.image[jjj, iii] < datamin:
                    datamin = self.image[jjj, iii]
                elif self.image[jjj, iii] > datamax:
                    datamax = self.image[jjj, iii]
                if self.image[jjj, iii] < ymin:
                    self.image[jjj, iii] = ymin
                elif self.image[jjj, iii] > ymax:
                    self.image[jjj, iii] = ymax

        print('Limit: Data Min: ', datamin, ', Data Max: ', datamax)
        print('Limited to  Min: ', ymin, ', Max: ', ymax)

        return

    def zero(self):
        """
        zeros() zeros the values for an image to make sure they are all finite
        """

        for iii in range(self.img_width):
            for jjj in range(self.img_height):
                self.image[jjj, iii] = 0
                self.weights[jjj, iii] = 0
        return

    def set_ij(self, iii, jjj, value, weight):
        """
        set_ij() sets a particular value of the grid
        """

        self.image[jjj, iii] = value
        self.weights[jjj, iii] = weight
        return

def main():
    """
    Define a main routine for testing of GridClass
    Only used if called without arguments
    """
    dpi = 2
    dpi = 1
    
    width = int(360)
    height = int(180)
    mywidth = int(width*dpi)
    myheight = int(height*dpi)
    xmin = 0.
    xmax = 360.
    ymin = -90.
    ymax = 90.
    FWHM = 13.0
    weight = 1.
    gridtype = 'GAL'
    mygrid = Grid(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, \
                      width=width, height=height, dpi=dpi, FWHM=FWHM,\
                      projection="-CAR", gridtype=gridtype)

    for iii in np.arange(xmin, xmax, 30):
        x = iii 
        for jjj in np.arange(ymin, ymax, 30):
            y = jjj
            z = y*x
            mygrid.convolve(x, y, z, weight)

    mygrid.normalize()
    mygrid.check()

    ticks = np.arange(0, mywidth, 30*dpi)
    x_ticks = xmin + (xmax-xmin)*ticks/mywidth
#y_ticks = ymin + (ymax-ymin)*ticks/myheight
    y_ticks = ymax - (ymax-ymin)*ticks/myheight
    plt.yticks(ticks, y_ticks)
    plt.xticks(ticks, x_ticks)

    plt.imshow(mygrid.image, interpolation='nearest')
    plt.title(gridtype + " Image")
    plt.show()

    for iii in np.arange(-mywidth, mywidth, 60*dpi):
        x = float(iii)
#        ix = mygrid.ii( x)
#        print "x -> ix: ", x, ix

if __name__ == "__main__":
    main()

