#Python Script to sort all observations and events into different directories
#HISTORY
#25Apr29 GIL Sort by file name, without reading the file
#25Feb25 GIL Add bigger directory range
#23Mar29 GIL look for all directories corresponding to a specific date
#23Feb15 GIL initial version
#
import sys
import os
import time
import datetime
import subprocess
import radioastronomy

dy = -1.

nargs = len( sys.argv)
if nargs < 2:
    print('sortdate: sort all of one days observations into UTC based dates')
    print('sortdate usage:')
    print('sortdate [-v] <date>')
    print('where')
    print('-v          optionally print verbose summary information')
    print('<date>      of observations to be sorted')
    print('            Ie 23Mar29')
    print('Will resort in all data and events being sorted for all telescopes')
    print(' Ie pi1-data-23Mar29, pi1-events-23Mar29, pi3-data-23Mar29 ...')
    print('  Only *.eve, *.ast, *.hot and *.cld files will be read. Others will be skipped')
    print("")
    print("Glen Langston, NSF -- 23 March 29")
    exit()

monthnames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", \
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ]
calendar = ""

dates = sys.argv[1:nargs]

# prepare to extract the prefix of the directories
prefix = ""

rs = radioastronomy.Spectrum()

# in case of first date is actually an argument
if dates[0] == '-v':
    verbose = True
    dates = dates[1:]
else:
    verbose = False

dates = sorted(dates)

nDates = len(dates)

dirtypes = [ "data", "events"]
count = 0

for date in dates:

    for ipi in range (1,13):
        for atype in dirtypes:
            
            dir = "pi%0d-%s-%s" % (ipi, atype, date)
            if verbose:
                print("dir %s " % (dir))
                
            if (os.path.exists(dir)): 
        
                files = os.listdir(dir)
                nfiles = len(files)
                print("%5d files in directory %s" % (nfiles, dir))
                sortcmd = "python3 ~/Research/analyze/sortday.py %s/*" % dir
                if verbose:
                    print("Sorting: %s" % (sortcmd))
                os.system(sortcmd)
                count = count + nfiles

if nDates > 1:
    print("%5d sortable files for date %s ... %s" % \
          (count, dates[0], dates[nDates-1]))
else:
    print("%5d sortable files for date %s" % (count, dates[0]))

