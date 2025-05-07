#Python Script to sort observations and events into different directories
#This version only uses the file name to determine the correct directory
#HISTORY  
#25Apr29 GIL initial version, based on sortdata.py
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
    print('sortday: sort spectra and events into different sub-directories')
    print('sortday usage:')
    print('sortday [-v] <filenames>')
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

if names[0].lower() == '-v':
    verbose = True
    names = names[1:]
else:
    verbose = False

count = 0
printcount = 0
movecount = 0
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

    # get the file name from the full path
    aname = parts[nparts-1]
    # current directory is adir
    adir = parts[0]
    if prefix == "":
        dirparts = adir.split("-")
        ndirparts = len(dirparts)
        if ndirparts != 3:
            print("Directory name must have 3 parts separated by '-'")
            print("Invalid Directory: %s" % adir)
            exit()
        prefix = "%s-%s-" % (dirparts[0], dirparts[1])
        inDir = adir
        if verbose:
            print("Prefix: %s; Current directory: %s" % (prefix,inDir))
        
    parts = aname.split('.')
    nparts = len(parts)
    if nparts == 1:
        printcount = printcount + 1
        print("%d: %s Only 1 part" % ( nparts, parts[0]))
        
    # if not a two part name with expected extension
    if nparts < 2:
        continue
# only summarize astronomy files
    if (parts[1] != 'ast') and (parts[1] != 'hot') and \
       (parts[1] != 'cld') and (parts[1] != 'eve'):
        continue
    if verbose and printcount < 10:
        printcount = printcount + 1
        print('File: ',filename)


    # expecting a name like "25-04-29T123456_789.eve"
    nameparts = parts[0].split('T')
    nnameparts = len(nameparts)
    if nnameparts == 1:
        print("Invalid file name: %s" % (parts[0]))
    # only sort if "T" is in the name
    if nnameparts != 2:
        continue
    # nameparts[0] now contains "25-04-29"
    dateparts = nameparts[0].split('-')
    ndateparts = len( dateparts)
    if ndateparts != 3:     # expecting YY MM DD, else skip
        continue
    imonth = int( dateparts[1]) - 1
    yy = dateparts[0]
    mm = dateparts[1]
    dd = dateparts[2]
        
    count = count + 1
    # check if the date changed
    if lastyy != yy or lastmm != mm or lastdd != dd:
        newdate = "%s%s%s" % (yy, monthnames[int(mm)-1], dd)
        lastyy = yy
        lastmm = mm
        lastdd = dd
        if calendar == "":
            calendar = newdate
            print("First Calendar Date: %s" % (calendar))
        newDir = "%s%s" % (prefix, newdate)

        mkdir = "mkdir %s 2> /dev/null" % (newDir)
        if not os.path.isfile( newDir):
            os.system(mkdir)
            if verbose:
                print("New Directory: %s " % (newDir))

    icount = int( count/100)
    if (newDir == inDir):
        
        if count < 3 or (100*icount) == count:
            print ("File %5d in the correct directory: %s" % (count, filename))
    else:
        if count < 3 or (100*icount) == count:
            print ("File %5d moved to new   directory: %s" % (count, newDir))
        mvDay = "mv %s %s 2> /dev/null" % (filename,newDir)
        movecount = movecount + 1
        if count < 5:
            print(mvDay)
        os.system(mvDay)
     
print("Moved %5d of %5d sortable files in  %s" % (movecount, count, inDir))

