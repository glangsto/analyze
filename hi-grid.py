# Create a velocity versus LST grid
# HISTORY
# 16NOV07 GIL initial version based on web example.

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma
from numpy.random import uniform, seed
# make up some randomly distributed data
seed(1234)
npts = 200
nLongitude = 360
lonMin = 0.
lonMax = 360.
nVelocity = 360
velMin = -300.
velMax = 300.
x = uniform(lonMin,lonMax,nLongitude)
y = uniform(velMin,velMax,nVelocity)
z = x*np.sin(x*y)
#print z
npts = len(z)
# define grid.
xi = np.linspace(lonMin,lonMax,nLongitude)
yi = np.linspace(velMin,velMax,nVelocity)
# grid the data.
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
print zi
# contour the gridded data, plotting dots at the randomly spaced data points.
CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
# plot data points.
plt.scatter(x,y,marker='o',c='b',s=5)
plt.xlim(lonMin,lonMax)
plt.ylim(velMin,velMax)
plt.title('griddata test (%d points)' % npts)
plt.show()
