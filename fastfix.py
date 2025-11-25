#Python Script to fix set of NSF events or spectra
#The list of header items to fix is given in help.
#These include El (elevation) and Az (azimuth)
#HISTORY
#25OCT10 GIL speed up processing
#25OCT08 GIL speed up processing
#
import sys
import radioastronomy
import numpy as np
import os

nargs = len( sys.argv)

# set magic values to identify new inputs
NOVALUE = -200.
newEl = NOVALUE
newAz = NOVALUE
newdEl = NOVALUE
newdAz = NOVALUE
newLat = NOVALUE
newLon = NOVALUE
newAlt = NOVALUE
newCen = NOVALUE
newBan = NOVALUE
newRefSample = NOVALUE
newRefChan = NOVALUE
newNChan = NOVALUE
newNTime = -200
newGain1 = NOVALUE
newGain2 = NOVALUE
newGain3 = NOVALUE
observer = ""
note = ""
telescope = ""
device = ""
# flag replacing file
replace = False
aFix = False

ifile = 1
iii = ifile

while iii < nargs:
    anarg = sys.argv[iii].upper()
    if str(anarg[0:3]) == "-EL":
        newEl = float( sys.argv[iii+1])
        iii = iii + 1
        print("New El: %7.2f" % (newEl))
        aFix = True
        ifile = ifile + 2
    if str(anarg[0:4]) == "-DEL":
        newdEl = float( sys.argv[iii+1])
        iii = iii + 1
        print("New El Offset: %7.2f" % (newdEl))
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-AZ":
        newAz = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Az: %7.2f" % (newAz))
        ifile = ifile + 2
        aFix = True
    if str(anarg[0:4]) == "-DAZ":
        newdAz = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Az Offset: %7.2f" % (newdAz))
        aFix = True
        ifile = ifile + 2
    if anarg[0:3] == "-LA":
        newLat = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Latitude : %13.9f" % (newLat))
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-LO":
        newLon = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Longitude: %13.9f" % (newLon))
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-AL":
        newAlt = float( sys.argv[iii+1])
        iii = iii + 1
        print("New Altitude : %6.3f" % (newAlt))
        ifile = ifile + 2
        aFix = True
    if anarg == "-NT":
        newNTime = np.int( sys.argv[iii+1])
        iii = iii + 1
        print("New Number of Samples: ", newNTime)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-OB":
        observer = sys.argv[iii+1]
        iii = iii + 1
        print("Observers: ", observer)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-NO":
        note = sys.argv[iii+1]
        iii = iii + 1
        print("Note: ", note)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-TE":
        telescope = sys.argv[iii+1]
        iii = iii + 1
        print("Telescope: ", telescope)
        ifile = ifile + 2
        aFix = True
    if anarg[0:4] == "-DEV":
        device = sys.argv[iii+1]
        iii = iii + 1
        print("Device: ", device)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-CE":
        newCen = float(sys.argv[iii+1])
        iii = iii + 1
        print("Center Frequency: ", newCen, " (MHz)")
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-BA":
        newBan = float(sys.argv[iii+1])
        iii = iii + 1
        print("Bandwidth ", newBan, " (MHz)")
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-RC":
        newRefChan = float(sys.argv[iii+1])
        iii = iii + 1
        print("Ref Chan  ", newRefChan)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-NC":
        newNChan = float(sys.argv[iii+1])
        iii = iii + 1
        print("N Channel ", newNChan)
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-RS":
        newRefSample = float(sys.argv[iii+1])
        iii = iii + 1
        print("Reference Sample ", newRefSample, " ")
        ifile = ifile + 2
        aFix = True
    if anarg[0:6] == "-GAIN1":
        newGain1 = float(sys.argv[iii+1])
        iii = iii + 1
        print("Gain 1: %7.2f" % (newGain1))
        ifile = ifile + 2
        aFix = True
    if anarg[0:6] == "-GAIN2":
        newGain2 = float(sys.argv[iii+1])
        iii = iii + 1
        print("Gain 2: %7.2f" % (newGain2))
        ifile = ifile + 2
        aFix = True
    if anarg[0:6] == "-GAIN3":
        newGain3 = float(sys.argv[iii+1])
        iii = iii + 1
        print("Gain 3: %7.2f" % (newGain3))
        ifile = ifile + 2
        aFix = True
    if anarg[0:3] == "-RE":
        replace = True 
        print("Replacing original file")
        ifile = ifile + 1
        aFix = True
    iii = iii + 1

#print "Afix: ",aFix, iii
# if nothing to fix, give help
if aFix == False:
    print("FASTFIX: Fix observing file parameters")
    print("Usage: Fix [-el elevation] [-az azimuth]... <file 1> [<file 2>] ... [<file N>]")
    print("Where optionally the following paramters may be fixed")
    print(" -az  Telescope azimuth in degrees")
    print(" -el  Telescope elevation in degrees")
    print(" -daz Telescope azimuth offset in degrees")
    print(" -del Telescope elevation offset in degrees")
    print(" -lat Telescope latitude in degrees")
    print(" -lon Telescope longitude in degrees")
    print(" -alt Telescope altitude above sea level in meters")
    print(" -cen Center Frequency (MHz)")
    print(" -ban Bandwidth (MHz)")
    print(" -rch Reference channel  (1 based)")
    print(" -nch Number of channels (1 based)")
    print(" -obs Observers names (ascii)")
    print(" -tel Telescope name (ascii)")
    print(" -not Note describing observation (ascii)")
    print(" -dev Software Defined Radio used for observation (ascii)")
    print(" -nt  Number of time samples in the observations")
    print(" -rs  Index of the Reference Time Sample")
    print(" -re  Replace original file with revised header")
    exit()

nFix = 0
nfiles = nargs-ifile
# create the default file structure
#rs = radioastronomy.Spectrum()

tempname = /tmp/fastfix.tmp
outname = ""

for iii in range(nfiles):

    filename = sys.argv[iii+ifile]

    if not os.path.isfile( filename):
        print("Input argument: %s is not a file:" % (filename))
        continue

    try:
        file_size = os.path.getsize(filename)
        if file_size < 2000:
            print(" Input file is too small (%d): %s" % (file_size, filename))
            continue
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
    except PermissionError:
        print(f"Error: Permission denied to access the file '{filename}'.")
    
    outfile = open(outname, 'w')
    infile = open(filename, 'f')

    if verbose:
        print filename

    try:
        with open(filename, 'r', encoding='ascii') as infile, \
             open(outname, 'w', encoding='ascii') as outfile:
            for line in infile:
                keyvalue = '# NOTEA     '
                if keyvalue in line:
                    if verbose:
                        
        outfile.write(outline)
        self.noteB = self.noteB.replace('\n', '')
        self.noteB = self.noteB.strip()
        keyvalue = '# NOTEB     '
        outfile.write(outline)
        self.observer = self.observer.replace('\n', '')
        self.observer = self.observer.strip()
        keyvalue = '# OBSERVER  '
        outfile.write(outline)
        self.device = self.device.replace('\n', '')
        self.device = self.device.strip()
        keyvalue = '# DEVICE    '
        outfile.write(outline)
        self.datadir = self.datadir.replace('\n', '')
        self.datadir = self.datadir.strip()
        keyvalue = '# DATADIR   '
        outfile.write(outline)
        self.site = self.site.replace('\n', '')
        self.site = self.site.strip()
        keyvalue = '# SITE      '
        outfile.write(outline)
        self.city = self.city.replace('\n', '')
        self.city = self.city.strip()
        keyvalue = '# CITY      '
        outfile.write(outline)
        self.region = self.region.replace('\n', '')
        self.region = self.region.strip()
        keyvalue = '# REGION    '
        outfile.write(outline)
        self.country = self.country.replace('\n', '')
        self.country = self.country.strip()
        keyvalue = '# COUNTRY   '
        outfile.write(outline)
        self.telType = self.telType.replace('\n', '')
        self.telType = self.telType.strip()
        keyvalue = '# TELTYPE   '
        outfile.write(outline)
        self.frame = self.frame.replace('\n', '')
        self.frame = self.frame.strip()
        keyvalue = '# FRAME     '
        outfile.write(outline)
        ngains = len(self.gains)
        if ngains > 0:
        keyvalueline = '# GAIN1     '
            outfile.write(outline)
        if ngains > 1:
        keyvalueline = '# GAIN2     '
            outfile.write(outline)
        if ngains > 2:
        keyvalueline = '# GAIN3     '
            outfile.write(outline)
        if ngains > 3:
        keyvalueline = '# GAIN4     '
            outfile.write(outline)
        keyvalue = '# Count     '
        outfile.write(outline)
        # match SETI/GUPPI KEYWORDS
        # https://www.cv.nrao.edu/~pdemores/GUPPI_Raw_Data_Format/
#       keyvaluee = '# CenterFreq'
        keyvalue = '# REFFREQ   '
        outfile.write(outline)
        keyvalue = '# OBSFREQ   '
        outfile.write(outline)
#       keyvaluee = '# Bandwidth '
        keyvalue = '# OBSBW     '
        outfile.write(outline)
        keyvalue = '# Duration  '
        outfile.write(outline)
        keyvalue = '# DeltaX    '
        outfile.write(outline)
        keyvalue = '# TSYS      '
        outfile.write(outline)
        keyvalue = '# TRX       '
        outfile.write(outline)
        keyvalue = '# TRMS      '
        outfile.write(outline)
        keyvalue = '# TINT      '
        outfile.write(outline)
        keyvalue = '# KPERC     '
        outfile.write(outline)
        keyvalue = '# GAINFACT  '
        outfile.write(outline)
        keyvalue = '# BUNIT     '
        outfile.write(outline)
        keyvalue = '# NCHAN     '
        outfile.write(outline)
        keyvalue = '# NSPEC     '
        outfile.write(outline)
        keyvalue = '# NTIME     '
        outfile.write(outline)
        keyvalue = '# NSAMPLES  '
        outfile.write(outline)
        keyvalue = '# EPEAK     '
        outfile.write(outline)
        keyvalue = '# ERMS      '
        outfile.write(outline)
        keyvalue = '# EMJD      '
        outfile.write(outline)
        keyvalue = '# REFCHAN   '
        outfile.write(outline)
        keyvalue = '# REFSAMPL  '
        outfile.write(outline)
        nave = self.nave
        keyvalue = '# NAVE      '
        outfile.write(outline)
        nmedian = self.nmedian
        keyvalue = '# NMEDIAN   '
        outfile.write(outline)
        keyvalue = '# Fft_rate  '
        outfile.write(outline)
        strnow = self.utc.isoformat()
        dates = strnow.split('T')
        datestr = dates[0] + ' ' + dates[1]
        keyvalue = '# UTC       '
        outfile.write(outline)
        keyvalue = '# SECONDS   '
        outfile.write(outline)
        lststr = angles.fmt_angle(self.lst/15., s1=":", s2=":", pre=3)  # convert to hours
        keyvalue = '# LST       '
        outfile.write(outline)
        keyvalue = '# AZ        '
        outfile.write(outline)
        keyvalue = '# EL        '
        outfile.write(outline)
        anglestr = angles.fmt_angle(float(self.tellon), s1=":", s2=":")
        keyvalue = '# TELLON    '
        outfile.write(outline)
        anglestr = angles.fmt_angle(float(self.tellat), s1=":", s2=":")
        keyvalue = '# TELLAT    '
        outfile.write(outline)
        keyvalue = '# TELALT    '
        outfile.write(outline)
        rastr = angles.fmt_angle(self.ra/15., s1=":", s2=":", pre=3) # convert to hours
        keyvalue = '# RA        '
        outfile.write(outline)
        decstr = angles.fmt_angle(self.dec, s1=":", s2=":")
        keyvalue = '# DEC       '
        outfile.write(outline)
        lonstr = angles.fmt_angle(self.gallon, s1=":", s2=":", pre=2)
        keyvalue = '# GALLON    '
        outfile.write(outline)
        latstr = angles.fmt_angle(self.gallat, s1=":", s2=":", pre=2)
        keyvalue = '# GALLAT    '
        outfile.write(outline)
        altstr = angles.fmt_angle(self.altsun, s1=":", s2=":", pre=1)
        keyvalue = '# ALT_SUN   '
        outfile.write(outline)
        az_str = angles.fmt_angle(self.az_sun, s1=":", s2=":", pre=1)
        keyvalue = '# AZ_SUN    '
        outfile.write(outline)
        keyvalue = '# ETAA      '
        outfile.write(outline)
        keyvalue = '# ETAB      '
        outfile.write(outline)
        keyvalue = '# POLANGLE  '
        outfile.write(outline)
        keyvalue = '# TELSIZEAM '
        outfile.write(outline)
        keyvalue = '# TELSIZEBM '
        outfile.write(outline)
        keyvalue = '# AST_VERS  '
        outfile.write(outline)
                
                if search_string in line:
                    modified_line = line.replace(search_string, replacement_string)
                    outfile.write(modified_line)
                else:
                    outfile.write(line)
        print(f"File '{input_filename}' processed. Output written to '{output_filename}'.")
    except FileNotFoundError:
        print(f"Error: The file '{input_filename}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    if newAz != NOVALUE:
        telaz = newAz
    if newEl != NOVALUE:
        telel = newEl
    if newdAz != NOVALUE:
        teldaz = newdAz
    if newdEl != NOVALUE:
        teldel = newdEl
    if newLat != NOVALUE:
        tellat = newLat
    if newLon != NOVALUE:
        tellon = newLon
    if newAlt != NOVALUE:
        telelev = newAlt
    if newGain1 != NOVALUE:
        gains[0] = newGain1
    if newGain2 != NOVALUE:
        gains[1] = newGain2
    if newGain3 != NOVALUE:
        gains[2] = newGain3
    if newAlt != NOVALUE:
        telelev = newAlt
    if observer != "":
        observer = observer
    if device != "":
        device = device
    if note != "":
        noteA = note
    if telescope != "":
        site = telescope
    if newNTime > 0:
        nTime = newNTime
    if newBan > 0:
        bandwidthHz = newBan * 1.E6
    if newCen > 0:
        centerFreqHz = newCen * 1.E6
    if newRefSample != NOVALUE:
        refSample = newRefSample
    if newRefChan != NOVALUE:
        refChan = newRefChan
    if newNChan != NOVALUE:
        nChan = newNChan

    parts = filename.split('/')
    nparts = len(parts)
    # get the file name without directory name
    aname = parts[nparts-1]
    filepart = aname
    # if no directory name
    if nparts == 1:
        dirname = "./"   # use current directory
    elif nparts == 2:
        dirname = parts[0]
    else:
        dirname = "./"
        for i in range(nparts-1):
            dirname = dirname + parts[i] + "/"
#    print "Directory: ", dirname

    extension = ""
    # if replacing original file
    if replace:
        try:
            os.remove(tempname)
        except:
            print("Cound not remove file: ",filename)
        outname = filename
        parts = filename.split('.')
        # last bit of file name is the extension
        nparts = len(parts)
        extension = parts[nparts-1]
    else:
        parts = filepart.split('.')
        nparts = len(parts)
        if (nparts == 2):
            extension = parts[1]
            outname = parts[0] + "-fix." + extension
        elif (nparts == 1):
            outname = parts[0] + "-fix"
        else:
            extension = parts[nparts-1]
            outname = parts[0] + "-fix." + extension
    
# now if a spectrum, fix name for elevation above zero 
# Spectra do not have time series.
    if rs.nTime <= 0:
        parts = outname.split('.')
        nparts = len(parts)
        # last part of file name is file type
        filetype = parts[nparts-1]
        # special case of a notes files
        if extension == "not":
            filetype = extension
        else:
            if rs.telel > 0.:
                filetype = 'ast'
            else:
                filetype = 'hot'
                
        if nparts == 2:
            outname = parts[0] + "." + filetype
        elif nparts == 1:
            outname = parts[0]
        else:  # else multiple name parts separated by "." 
            outname = ""
            for iii in range(nparts-1):
                outname = outname + parts[iii] + "."
            outname = outname + filetype
#    print "Output file name: ", outname
    rs.write_ascii_file( dirname, outname)
    nFix = nFix + 1

print("Fixed %d Files" % (nFix))

def main():
    """
    Test file if no arguments
    """

    # Example usage:
    input_file = "input.txt"
    output_file = "output.txt"
    string_to_find = "old_text"
    string_to_replace_with = "new_text"

    # Create a dummy input file for demonstration
    with open(input_file, 'w', encoding='ascii') as f:
        f.write("This is a line with old_text.\n")
        f.write("Another line without the string.\n")
        f.write("A third line containing old_text multiple times: old_text, old_text.\n")

        modify_file_line(input_file, output_file, string_to_find, string_to_replace_with)

    return

