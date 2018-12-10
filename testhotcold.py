#Python Script to plot raw NSF record data.
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#18MAR05 GIL implement folding option
#18FEB06 GIL keep track of first az,el
#17AUG18 GIL plot the last spectrum in the series
#17AUG17 GIL Note elevation range
#16Oct07 GIL check for changes in frequency and/or gain
#16AUG29 GIL make more efficient
#16AUG16 GIL use new radiospectrum class
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime
import statistics
import radioastronomy
import copy
import interpolate
import radioastronomy

# default values
avetimesec = 3600.
# some SDRs put spike in center of spectrum; indicate spike flagging here
flagCenter = False
flagCenter = True
# put list of RFI features here, for interpolation later
flagRfi = False
flagRfi = True
# put your list of known RFI features here.  Must have at least two.
linelist = [1400.00, 1420.0]  # RFI lines in MHz
linewidth = [5, 5]

nargs = len(sys.argv)
#first argument is the averaging time in seconds
timearg = 1
namearg = 2

scalefactor = 1e8
xallmax = -9.e9
xallmin = 9.e9
yallmax = -9.e9
yallmin = 9.e9
# velocities for fitting baselines
minvel = -450.
minvel = -280.
maxvel = 210.
# currently used
maxvel = 160.
minvel = -160.
# currently used
maxvel = 200.
minvel = -180.
# currently used
maxvel = 180.
minvel = -550.

c = 299792.458  # (v km/sec)
nuh1 = 1420.4057517667 # neutral hydrogen frequency (MHz)
thot = 285.0  # define hot and cold
#thot = 272.0  # 30 Farenheit = 272 K
tcold = 10.0
tmin = 20.0 
tmax = 999.0 # define reasoanable value limits

# prepare to remove a linear baseline

xa = 200
#xb = 1024-xa
xb = 1024-200

#for symbol, value in locals().items():
#    print symbol, value

dt = datetime.timedelta(seconds=0.)

#first argument is the averaging time in seconds
names = sys.argv[namearg:]
names = sorted(names)
nFiles = len(names)
rs = radioastronomy.Spectrum()

filename = names[0]

rs.read_spec_ast(filename)

vels = [-200., 100]
print "vels : ", vels
freqs = rs.vel2freq( vels, nuh1*1.E6)

print "freqs: ", freqs

chans = rs.freq2chan( freqs)
print "chans: ", chans


chans = rs.vel2chan( vels, nuh1*1.E6)
print "chans: ", chans

freqs = rs.chan2freq( chans)
print 'freqs: ', freqs

vels = rs.chan2vel( chans, nuh1*1.E6)
print 'vels : ', vels
