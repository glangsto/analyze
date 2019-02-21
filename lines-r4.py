# script to estimate linear molecule line strenght as a function of temperature
# Glen Langston
# HISTORY
# 17OCT02 GIL Complete full line frequency calculation
# 17OCT01 GIL Initial version
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma

# define important constants
BMHz = 564.0011         #Rotational constant MHz
sigma = 2.              #energy separation factor for frequency
Bcm = 29997.1/BMHz      #convert to wavelength cm.    
Bcm1 = 1./Bcm           #convert to inverse cm
kb300Koverhc = 209.7    # cm-1 for 300 k temp
T = 1.                 # physical temperature
kboverhc = kb300Koverhc * T / 300. 
qrot = kboverhc / ( Bcm1 * sigma)
muA = 4.82
npts = 200

print 'KB/hc: ', kboverhc
print 'q_rot: ', qrot
xx = np.zeros((npts))
nu = np.zeros((npts))
yy = np.zeros((npts))
ee = np.zeros((npts))
################
# From CDMS: 
#Intensity Calculations
#E" and E' the lower and upper state energy, respectively, and 
#Qrs the rotation- spin partition function at the temperature T

NtoEnergy = 4.16231E-5

maxy = 0.
ee[0] = np.exp( -0.)
yy[0] = 0.
ysum = yy[0]
for iii in range(1,npts):
    xx[iii] = float(iii)
    nu[iii] = xx[iii] * BMHz * sigma
    ee[iii] = np.exp( - xx[iii] / qrot)
    sg = ((2.*xx[iii]) + 1.)
    yy[iii] = NtoEnergy * sg * (ee[iii-1] - ee[iii])
    ysum = ysum + yy[iii]
    if yy[iii] > maxy:
        maxy = yy[iii]
print 'Computed partition function for ', npts, ' transitions'
print 'Sum : ', ysum
yy = yy / ysum

mytitle='Test of Plotting Partition Functions'
fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
plt.plot(xx, yy, 'b', linestyle='-', label='HC7N', lw=4)
plt.title(mytitle, fontsize=16)
plt.xlabel('Frequency (MHz)', fontsize=16)
plt.ylabel('Line Intensity (Relative)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
