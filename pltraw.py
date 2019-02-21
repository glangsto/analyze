#Python Script to plot raw NSF record data.
#
#import matplotlib.pyplot as plt
#plot the raw data from the observation
#HISTORY
#15AUG30 add option to plot range fo values
#15JUL01 GIL Initial version
#
import pylab
import numpy as np
import sys
import statistics

dy = -1.
scalefactor=1e10

print 'Argument List:', str(sys.argv)
nargs = len( sys.argv)

linestyles = ['-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.']
colors = ['-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g','-b','-r','-g']

pylab.figure()
pylab.hold(True)

xallmax = -9.e9
xallmin =  9.e9
yallmax = -9.e9
yallmin =  9.e9
az = '-1.'
el = '-90'
glat = 0
glon = 0
count = '0'
note = ''

for iii in range(1, min(nargs,10)):

    filename = sys.argv[iii]
#    print(iii,filename)
# Read the file. 
    f2 = open(filename, 'r')
# read the whole file into a single variable, which is a list of every row of the file.
    lines = f2.readlines()
    f2.close()

# initialize some variable to be lists:
    x1 = []
    y1 = []

    datacount = 0
    linecount = 0
# scan the rows of the file stored in lines, and put the values into some variables:
    for line in lines:
        parts = line.split()
#        if linecount == 2:
#            print('%s' % (line))
        linecount = linecount + 1
        if linecount == 3:
            note = line

# if a comment or parameter lien
        if line[0] == '#':
            if len(parts) < 4:
                continue
            if parts[1] == 'AZ':
                az = parts[3]
            if parts[1] == 'EL':
                el = parts[3]
            if parts[1] == 'LON':
                lon = parts[3]
                lonparts = lon.split(':')
                lonfloat = float(lonparts[0]) + float(lonparts[1])/60. + float(lonparts[2])/3600.
# if latitude parse and turn into a float
            if parts[1] == 'LAT':
                lat = parts[3]
                latparts = lat.split(':')
                if float(latparts[0] > 0):
                    latfloat = float(latparts[0]) + float(latparts[1])/60. + float(latparts[2])/3600.
                else:
                    latfloat = -float(latparts[0]) + float(latparts[1])/60. + float(latparts[2])/3600.
                    latfloat = - latfloat
            continue
        datacount = datacount+1
        p = line.split()
#        print(p)
        x1.append(float(p[1])*1.E-6)
        y1.append(float(p[2])*scalefactor)
    
    xv = np.array(x1)
    yv = np.array(y1)
    
    label = '%s, AZ,EL: %5s,%5s, Lon,Lat=%5.1f,%5.1f' % (filename,az,el,lonfloat,latfloat)
    xmin = min(xv)
    xmax = max(xv)
    xallmin = min(xmin,xallmin)
    xallmax = max(xmax,xallmax)
    ymin = min(yv)
    ymax = max(yv)
    ymed = statistics.median(yv)

#    ymed = (yv[300] + yv[700])*.5
    print('File: %s' % filename) 
    print(' Max: %9.1f  Median: %9.1f SNR: %6.2f ; %s %s' % (ymax, ymed, ymax/ymed, count, label)) 
    print('%s' % note)
    yallmin = min(ymin,yallmin)
    yallmax = max(ymax,yallmax)
    pylab.xlim(xallmin,xallmax)
#    pylab.ylim(1.1*yallmin,.9*yallmax)
#    pylab.ylim(.2*yallmax,1.2*yallmax)
    pylab.ylim(0,1.5*yallmax)
#    pylab.ylim(yallmin,150)
#    pylab.ylim(500, 2000)
    pylab.title(note)
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Intensity (Counts)')
#    print( 'Plotting ',datacount,' channels from ',filename)
# now, plot the data:
    pylab.plot(xv, yv, colors[iii-1], linestyle=linestyles[iii-1],label=label)
    pylab.legend(loc='upper right')
pylab.show()
