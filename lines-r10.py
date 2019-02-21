# script to estimate linear molecule line strenght as a function of temperature
# Glen Langston
# HISTORY
# 17OCT05 GIL output in a form similar to CDMS
# 17OCT02 GIL merge different notations
# 17OCT02 GIL Complete full line frequency calculation
# 17OCT01 GIL Initial version
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma

# Physica Temperature of this simulation in K
molecule = 'HC7N'
T = 9.375               # Kelvins
T = 4.7
#T = 3.                  # Kelvins
#T=50.

#constants for molecule (HC7N)
BMHz = 564.00112        #Rotational constant Frequency (MHz) from McCarthy
BMHz = 564.0011         #Rotational constant Frequency (MHz) from McCarthy
BMHz = 564.00113        #Rotational constant Frequency (MHz) from McCarthy
DMHz = 4.02             #stretch rotational constant (Hz)
eqQMHz = -4.29          # stretch reduction in frequency 
sigma = 2.              #energy separation factor for frequency
muA = 4.82              #A axis dipole moment

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
Dcm = DMHz/c               #stretch frequency in wave numbers
mua = 4.82

print 'k_b/hc: ', kboverhc, ' 1/sK'
Jmost = np.sqrt( kb * T / (2. * h * c * Bcm)) - 0.5
print 'For %s, at Temperature %7.2f K, the most populated  quantum state is J = %5.1f' % ( molecule, T, Jmost)
# precompute
nfactor = np.exp(-h*2.*BMHz*1.E6 / (kb*T))
print 'Between any to J levels, the exponential factor of level population is: %8.4e' % ( nfactor)
print 'Lowest frequency corresponds to temperature: %9.4f (K)' % ( 2.*BMHz*1.E6*h/kb)

xx = np.zeros((npts))
nu = np.zeros((npts))
Ej = np.zeros((npts))
Nj = np.zeros((npts))
Bj = np.zeros((npts))
Bnu = np.zeros((npts))
Qj = np.zeros((npts))
Jnu = np.zeros((npts))

################
# From CDMS: 
#Intensity Calculations
#E" and E' the lower and upper state energy, respectively, and 
#Qrs the rotation- spin partition function at the temperature T

NtoEnergy = 4.16231E-5

Ej[0] = 0
Nj[0] = 1.
nu[0] = 2. * BMHz * 1.E6
Bj[0] = 0.
Qj[0] = 1.
qsum = Qj[0]
bsum = Bj[0]
jmax = 0
Bmax = 0.
Jnu[0] = 0
mname = '%s_%.1fK' % (molecule, T)
for jj in range(0,npts):
    xx[jj] = float(jj)
    stretch = DMHz * xx[jj] * xx[jj] * (xx[jj] + 1.) * (xx[jj] + 1.)
    stretch = 4. * DMHz * (xx[jj]+1.) * (xx[jj] + 1.) * (xx[jj] + 1.)
    # note the transition energy between two levels and the energy of the level itself are different
    nu[jj] = (2.* xx[jj]*BMHz*1.E6 ) - stretch  # in Hz
    Ej[jj] = h * ( (xx[jj] * (xx[jj] + 1.) * BMHz * 1.E6) - stretch)             # Energy g cm^2/sec^2
    gj = ((2. * xx[jj]) + 1.)
    Nj[jj] = gj * np.exp( - Ej[jj]/(kb*T))
    Qj[jj] = gj * np.exp( - Ej[jj]/(kb*T))
    denom = (np.exp (h * Nj[jj] / (kb * T)) - 1.)
    if denom > 1.e-25:
        Jnu[jj] = (nu[jj] * h / kb) / denom
    else: 
        Jnu[jj] = 0.0
    denom = (np.exp(h*nu[jj]/(kb*T)) - 1.)
    if denom > 1e-28:
        Bj[jj] = (2. * h * nu[jj]*Nj[jj]*nu[jj]/(c*c))/denom
    else:
        Bj[jj] = 0.
    qsum = qsum + Qj[jj]
    bsum = bsum + Bj[jj]

    if Bj[jj] > Bmax:
        Bmax = Bj[jj]
        jmax = jj

for jj in range( 1, npts):
    if jj < 50:
#        print '%3d: Nu: %11.4f, N: %9.3e, Line Strength: %10.3e' % ( jj, nu[jj]*1.E-6, Nj[jj], Bj[jj])
        amp = np.log( Bj[jj]/bsum)
        print '%12.4f %12.4f %12.4f 3 %10.6f 123456 12345 67 89 %3d %3d %s' % \
            (nu[jj]*1.E-6, 0.001, amp, T, jj, jj-1, mname)

Bj = 100. * Bj / bsum
Bmax = 100. * Bmax / bsum
print 'Computed partition function for %3d transitions: %15.5f' % (npts, qsum)
print 'Maximum model intensity is for transition %3d' % (jmax)

mytitle='Molecular Line Relative Strengths for %8.3f K' % (T)
fig, ax1 = plt.subplots(figsize=(10, 6))
plt.hold(True)
# convert back to MHz for plotting
nu = nu * 1.E-6
#plt.plot(nu, Bj, 'b', linestyle='-', label=molecule, lw=4)
plt.ylim( 0., Bmax*1.1)
for jj in range( 1, npts):
    # search from high frequency end
    kk = npts - jj
    # if line intensity is signficiant
    if Bj[kk] > 0.001:
        # use this as the max x
        maxnu = nu[kk]
        break

plt.xlim( 0., maxnu)
# now plot each transition
plt.plot( [nu[1], nu[1]], [0.,Bj[1]], 'k-', lw=2, label=molecule)
for jj in range( 2, npts):
    plt.plot( [nu[jj], nu[jj]], [0.,Bj[jj]], 'k-', lw=2)

# finally plot survey frequency range:
#plt.axvline(x=38.9E3, color='r', linestyle='-.', lw=2)
#plt.axvline(x=48.0E3, color='r', linestyle='-.', lw=2)
# Set the tick labels font
for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
    label.set_fontsize(14)
#    label.set_fontname('Arial')
plt.axvspan(38.9E3, 48.0E3, color='b', alpha=0.1)
plt.title(mytitle, fontsize=16)
plt.xlabel('Frequency (MHz)', fontsize=16)
plt.ylabel('Line Intensity (% of Total Emission)', fontsize=16)
plt.legend(loc='upper right')
plt.show()
