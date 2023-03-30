#Python Script to sort observations and events into different directories
#HISTORY
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
    print('sortdata: sort spectra and events into different sub-directories')
    print('sortdata usage:')
    print('sortdata [-v] <filenames>')
    print('where')
    print('-v          optionally print verbose summary information')
    print('<filenames> list of file names to sort')
    print('  Only *.eve, *.ast, *.hot and *.cld files will be read. Others will be skipped')
    print("")
    print("Glen Langston, NSF -- 23 Feb 15")
    exit()

monthnames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", \
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ]
calendar = ""
batcmd="/bin/date +%Z"
timezone = subprocess.check_output(batcmd, shell=True)
parts = timezone.split()
# fix strange problem with zone sometimes coming back in quotes like b'EST'
try:
    timezone = str(parts[0], 'UTF-8')
except:
    timezone = str(parts[0])

# prepart to compute local time
#utcOffset = datetime.datetime.utcnow() - datetime.datetime.now()
#utcOffsetSecs = utcOffset.total_seconds() + 1.
#utcOffsetSecs = int(utcOffsetSecs)
#utcOffsetHours = utcOffsetSecs/3600.
#print("Time Zone: %s" % (timezone))

nsame = 0
nread = 0

names = sys.argv[1:nargs]
names = sorted(names)

# prepare to extract the prefix of the directories
prefix = ""

rs = radioastronomy.Spectrum()

if names[0] == '-v':
    verbose = True
    names = names[1:]
else:
    verbose = False

count = 0
lastyy = ""
lastmm = ""
lastdd = ""

for filename in names:

#    print filename
    parts = filename.split('/')
    nparts = len(parts)
    if nparts < 2:
        print("Names must include a directory, to extract prefix")
        print("Invalid: %s" % filename)
        exit()
    
    aname = parts[nparts-1]
    adir = parts[0]
    if prefix == "":
        dirparts = adir.split("-")
        ndirparts = len(dirparts)
        if ndirparts != 3:
            print("Directory name must have 3 parts separated by '-'")
            print("Invalid Directory: %s" % adir)
            exit()
        prefix = "%s-%s-" % (dirparts[0], dirparts[1])
        if verbose:
            print("Prefix: %s" % (prefix))
        inDir = adir
        
    parts = aname.split('.')
    nparts = len(parts)
    # if not a two part name with expected extension
    if nparts < 2:
        continue
# only summarize astronomy files
    if (parts[1] != 'ast') and (parts[1] != 'hot') and \
       (parts[1] != 'cld') and (parts[1] != 'eve'):
        continue
    if verbose:
        print('File: ',filename)
    rs.read_spec_ast( filename)
    if calendar == "":
        utcstr = str(rs.utc)
        strparts = utcstr.split()
        calendar = strparts[0]
        # print("First Calendar Date: %s" % (calendar))
        
    count = count + 1
    adatetime = str(rs.utc)
    parts = adatetime.split(' ')
    date  = parts[0]
    time  = parts[1]
    parts  = time.split('.')
    time = parts[0]

    dateparts = date.split("-")
    yy = dateparts[0]
    yy = yy[2:]
    mm = dateparts[1]
    dd = dateparts[2]
    # check if the date changed
    if lastyy != yy or lastmm != mm or lastdd != dd:
        newdate = "%s%s%s" % (yy, monthnames[int(mm)-1], dd)
        lastyy = yy
        lastmm = mm
        lastdd = dd
        newDir = "%s%s" % (prefix, newdate)

        mkdir = "mkdir %s 2> /dev/null" % (newDir)
        if not os.path.isfile( newDir):
            os.system(mkdir)
            if verbose:
                print("New Directory: %s " % (newDir))

    mvData = "mv %s %s 2> /dev/null" % (filename,newDir)
    icount = int( count/100)
    if (newDir == inDir):
        
        if count < 3 or (100*icount) == count:
            print ("File %5d in the correct directory: %s" % (count, filename))
    else:
        if count < 3 or (100*icount) == count:
            print ("File %5d moved to new   directory: %s" % (count, newDir))
        os.system(mvData)
     
print("%5d sortable files in  %s" % (count, inDir))

