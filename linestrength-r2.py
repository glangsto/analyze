# script to estimate linear molecule line strenght as a function of temperature
# Glen Langston
# HISTORY
# 17OCT01 GIL Initial versoin
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma

# define important constants
npts = 300
NtoEnergy = .2

x = uniform(1.,100.,npts)
xx = np.zeros((npts))
yy = np.zeros((npts))
y = ((2.*x)+1.)*np.exp(-(NtoEnergy*x))
z = ((2.*x)+1.)*np.exp(-(NtoEnergy*x))
zsum = 0.
ysum = 0.
for iii in range(0,npts):
    zsum = zsum + z[iii]
    xx[iii] = float(iii)
    ee = NtoEnergy * xx[iii]
    yy[iii] = ((2.*xx[iii]) + 1.) * np.exp(-(ee))
    ysum = ysum + yy[iii]
print 'Computed partition function for ', npts, ' transitions'
print 'Sum : ', ysum
yy = yy / ysum
y = z/zsum
maxy = np.max(y)
# define grid.
xi = np.linspace(1.,100.,100)
yi = np.linspace(0.,1.1*maxy,100)
# grid the data.
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
# contour the gridded data, plotting dots at the randomly spaced data points.
CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
# plot data points.
plt.scatter(x,y,marker='o',c='b',s=5)
plt.xlim(0.,100.)
plt.ylim(0.,1.1*maxy)
plt.title('Line Strength Calculation test (%d points)' % npts)
plt.show()

mytitle='Test of Plotting Partition Functions'
fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
plt.plot(xx, yy, 'b', linestyle='-', label='HC7N', lw=4)
plt.title(mytitle, fontsize=16)
plt.xlabel('Frequency (MHz)', fontsize=16)
plt.ylabel('Line Intensity (Relative)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
