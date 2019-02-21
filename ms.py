#Python Script to median a series of spectra as well as 
#record the min and max in this series
#HISTORY
#18FEB16 GIL only output if nfiles > 1
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import statistics
import radioastronomy
import copy
import interpolate

dy = -1.

nargs = len(sys.argv)

linestyles = ['-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--','-','-','-', '-.','-.','-.','--','--','--']
colors =  ['-b','-r','-g', '-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']
nmax = len(colors)
scalefactor = 1e8
xallmax = -9.e9
xallmin = 9.e9
yallmax = -9.e9
yallmin = 9.e9

c = 299792.  # (v km/sec)
nuh1 = 1420.4056 # neutral hydrogen frequency (MHz)
thot = 285.0  # define hot and cold
#thot = 272.0  # 30 Farenheit = 272 K
tcold = 10.0
tmin = 20.0 
tmax = 999.0 # define reasoanable value limits

#for symbol, value in locals().items():
#    print symbol, value

nplot = 0
nhot = 0
ncold = 0
lastfreq = 0.
lastbw = 0.
lastgain = 0.
lastel = 0.
lastaz = 0.
firstdate = ""
lastdate = ""
minel = 200.
maxel = -200.

# first read through all data and find hot load
names = sys.argv[1:]
names = sorted(names)
my_rs = []
count = 0
imin = count
imax = count

def median4( y0, y1, y2, y3):
    """
    median 4 arrays
    """
    n = len(y1)
    yout = np.zeros(n)
    y4 = np.zeros(4)
    for ii in range(n):
        y4[0] = y0[ii]
        y4[1] = y1[ii]
        y4[2] = y2[ii]
        y4[3] = y3[ii]
        yout[ii] = statistics.median( y4)

    return yout

nnames = len( names)
meds = np.zeros(nnames)
for filename in names:

    rs = radioastronomy.Spectrum()
    rs.read_spec_ast(filename)
    rs.azel2radec()    # compute ra,dec from az,el
    yv = rs.ydataA
    n = len(yv)
    n4 = round(n/4)
    n8 = round(n/8)
    na = int(n4)
    nb = int(na+n8)
    nd = int(3*n4)
    nc = int(nd-n8)
    ymedab = statistics.median(yv[na:nb])    
    ymedcd = statistics.median(yv[nc:nd])    
    # finally average values on opposite sides of spectrum
    ymed = (ymedab + ymedcd) *.5
    # now record signal away from center and lines; related to Tsys
    meds[count] = ymed

    # next keep track of spectra min and max values
    if meds[imin] > ymed:
        imin = count
    if meds[imax] < ymed:
        imax = count
        
    my_rs.append(rs)
    count = count + 1

ymean = statistics.mean( meds) * scalefactor
ystd  = statistics.stdev( meds) * scalefactor
ymedian = statistics.median( meds) * scalefactor
print 'Spectra power median = %.2f, mean = %.2f +/- %.2f : ' % ( ymedian, ymean, ystd)

print 'Read ',count,' Files '
print 'Minimum spectrum is : ', names[imin], meds[imin] * scalefactor
print 'Maximum spectrum is : ', names[imax], meds[imax] * scalefactor

nmedians = int(nnames/4)

my_copy = my_rs

for jj in range(nmedians):

    my_meds = []
    nloop = len(my_copy)
    medcount = 0
    for ii in range(0,nloop-3,4):
        print 'II: ',ii
        rsa = my_copy[ii]
        print rsa.utc
        rsb = my_copy[ii+1]
        rsc = my_copy[ii+2]
        rsd = my_copy[ii+3]
        
        ya = rsa.ydataA
        yb = rsb.ydataA
        yc = rsc.ydataA
        yd = rsd.ydataA
        
        yout = median4( ya, yb, yc, yd)
        
        rsb.ydataA = yout
        
        my_meds.append( rsb)
        medcount = medcount + 1

    my_copy = my_meds
    if medcount < 1:
        break
    print 'Median Count: ', medcount
    
fig,ax1 = plt.subplots(figsize=(10,6))
fig.canvas.set_window_title("Median of Spectra")

rsmax = my_rs[imax]
ymax = rsmax.ydataA*scalefactor
rsmin = my_rs[imin]
ymin = rsmin.ydataA*scalefactor

ymedian = rsb.ydataA*scalefactor
xv = rsmax.xdata  * 1.E-6 # convert to MHz

rsb.write_ascii_file( "./", "Median.dat")

plt.plot(xv, ymax, colors[1], linestyle=linestyles[1], label='Max', lw=2)
plt.plot(xv, ymedian, colors[2], linestyle=linestyles[1], label='Median', lw=2)
plt.plot(xv, ymin, colors[3], linestyle=linestyles[1], label='Min', lw=2)
plt.xlabel('Frequency (MHz)',fontsize=16)
plt.ylabel('Intensity (Counts)', fontsize=16)
plt.legend(loc='upper right')
plt.show()


