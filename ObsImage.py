# test of plotting Horn observations on an all Sky map
# HISTORY
# 16NOV08 GIL initial version

import numpy
from PIL import Image
from pylab import *

# map in galactic coordinates
filename = "/Users/glangsto/Desktop/Research/milkyway_21cm.jpg"
# map in ra-dec coordinates
filename = "/Users/glangsto/Desktop/Research/HEC_visible_sky_location_HR_crop.png"
filename = "/Users/glangsto/Desktop/Research/HEC_visible_sky_location_crop.png"
filename = "/Users/glangsto/Desktop/Research/1000Mess.png"

# read image to array
fig, ax1 = plt.subplots(1,1)
plt.tight_layout()
im = array(Image.open(filename))
nx = im.shape[1]
ny = im.shape[0]
print im.shape, im.dtype
xticks = [0, nx/4, nx/2, 3*nx/4, nx-1]
yticks = [0, ny/4, ny/2, 3*ny/4, ny-1]
ax1.set_yticks(yticks)
ax1.set_xticks(xticks)
xlabels = ax1.set_xticklabels(('24h','18h','12h','6h','0h'))
ylabels = ax1.set_yticklabels(('90d','45d','0d','-45d','-90d'))
ax1.set_xlabel("Local Siderial Time (hours)")
ax1.set_ylabel("Declination (degrees)")

#, extent=[24,0,-90,90])
imshow(im)

# some points
x = [100,100,400,400]
y = [200,500,200,500]

# plot the points with red star-markers
#plot(x,y,'r*')

# line plot connecting the first two points
#plot(x[:2],y[:2])

# add title and show the plot
title('Plotting: %s' % (filename))
show()

