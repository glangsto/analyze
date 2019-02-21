# script to estimate linear molecule line strenght as a function of temperature
# Glen Langston
# HISTORY
# 17OCT02 GIL Complete full line frequency calculation
# 17OCT01 GIL Initial version
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma

# Physica Temperature of this simulation in K
molecule = 'HC7N'
T = 9.375               # Kelvins

#constants for molecule (HC7N)
BMHz = 564.00112        #Rotational constant Frequency (MHz) from McCarthy
D = 4.02                #stretch rotational constant (Hz)
eqQMHz = -4.29          # stretch reduction in frequency 
sigma = 2.              #energy separation factor for frequency
muA = 4.82

# define important constants
pi = 3.1415926          #pi
c = 2.99792458E10       #speed of light, cm/sec
hbar = 1.0456E-27       #Planck's constant /2 pi, cgs cm^2 g/sec
h = hbar*(2.*pi)        #Planck's constant 
kb =  1.3807E-16        #Boltzmann's constant cm^2 g / (s^2 K)
kboverhc = kb/(c*h)     

npts = 100

# some useful deduced constants
Bcm = BMHz*1.E6/c       #Frequency in wave numbers
Dcm = D/c               #stretch frequency in wave numbers
qrot = kboverhc / ( Bcm * sigma)

print 'k_b/hc: ', kboverhc
print 'q_rot: ', qrot
Jmost = np.sqrt( kb * T / (2. * h * c * Bcm)) - 0.5
print 'For %s, at Temperature %7.2f K, the most populated  quantum state is J = %5.1f' % ( molecule, T, Jmost)
# precompute
nfactor = np.exp(-h*BMHz*1.E6 / (kb*T))
print 'Between any to J levels, the exponential factor of level population is: %8.4e' % ( nfactor)
print 'Lowest frequency corresponds to temperature: %9.4f (K)' % ( 2.*BMHz*1.E6*h/kb)

xx = np.zeros((npts))
nu = np.zeros((npts))
yy = np.zeros((npts))
Ej = np.zeros((npts))
Nj = np.zeros((npts))
Bj = np.zeros((npts))

################
# From CDMS: 
#Intensity Calculations
#E" and E' the lower and upper state energy, respectively, and 
#Qrs the rotation- spin partition function at the temperature T

NtoEnergy = 4.16231E-5

maxy = 0.
Ej[0] = 0
Nj[0] = 1.
yy[0] = 0.
nu[0] = 2. * BMHz * 1.E6
Bj[0] = 0.
ysum = yy[0]
for jj in range(1,npts):
    xx[jj] = float(jj)
    stretch = 4. * Dcm * (xx[jj] + 1.) * (xx[jj] + 1.) * (xx[jj] + 1.)
    nu[jj] = (2.* xx[jj]*BMHz*1.E6 ) - stretch
    Ej[jj] = h * c * Bcm * xx[jj] * (xx[jj]+1.)
    gj = ((2. * xx[jj]) + 1.)
    Nj[jj] = gj * np.exp( - Ej[jj]/(kb*T))
    Bj[jj] = (2. * h * nu[jj]*Nj[jj]*nu[jj]/(c*c))/(np.exp(h*nu[jj]/(kb*T)) - 1.)
    yy[jj] = NtoEnergy * gj * (Nj[jj-1] - Nj[jj])
    ysum = ysum + yy[jj]
    if jj < 15:
        print '%3d: Nu: %11.4f, N: %9.3e, Line Strength: %10.3e' % ( jj, nu[jj]*1.E-6, Nj[jj], yy[jj])
    if yy[jj] > maxy:
        maxy = yy[jj]
                         
print 'Computed partition function for ', npts, ' transitions'
print 'Sum : ', ysum
yy = yy / ysum

mytitle='Test of Plotting Partition Functions'
fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
plt.plot(xx, yy, 'b', linestyle='-', label=molecule, lw=4)
plt.title(mytitle, fontsize=16)
plt.xlabel('Frequency (MHz)', fontsize=16)
plt.ylabel('Line Intensity (Relative)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
