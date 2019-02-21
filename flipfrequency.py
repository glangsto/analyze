"""
Class defining a Radio Frequency Spectrum
Includes reading and writing ascii files
"""

##################################################
# Imports
##################################################
import datetime
import numpy
import angles
import radioastronomy

class Spectrum(object):
    """
    Define a Radio Spectrum class for processing, reading and
    writing astronomical data.
    """
    def __init__(self):
        """
        initialize all spectrum class values
        many will be overwritten laters
        """
        noteA = ""
        noteB = ""
        gains = [1., 1., 1.]
        utc = datetime.datetime.utcnow()
        telType = "Pyramid Horn"
        refChan = MAXCHAN/2
        observer = "Glen Langston"
        site = "Moumau House"
        xdata = numpy.zeros(MAXCHAN)
        ydataA = numpy.zeros(MAXCHAN)
        ydataB = numpy.zeros(MAXCHAN)
        #now fill out the spectrum structure.
        self.writecount = 0
        self.count = int(0)          # count of spectra summed
        self.noteA = str(noteA)      # observing note A
        self.noteB = str(noteB)      # observing note B
        self.observer = str(observer)   # name of the observer
        self.site = str(site)       # name of the observing site
        self.gains = gains           # one or more gain parameters
        self.telaz = 0.              # telescope azimuth (degrees)
        self.telel = 0.    # telescope elevation (degrees)
        self.tellon = 0.   # geographic longitude negative = West (degrees)
        self.tellat = 0.   # geopgraphic latitude (degrees)
        self.telelev = 0.  # geographic elevation above sea-level (meteres)
        self.centerFreqHz = 1.0   # centerfrequency of the observation (Hz)
        self.bandwidthHz = 1.0   # sampleRate of the observation (Hz)
        self.deltaFreq = 1.0   # frequency interval between channels
        self.utc = utc   # average observation time (datetime class)
        self.lst = 0.    # local sideral time degrees, ie 12h = 180deg
        self.durationSec = 0.    # integrated observing time (seconds)
        self.telType = str(telType) # "Horn, Parabola  Yagi, Sphere"
        # define size of horn or antenna (for parabola usuall A = B)
        self.telSizeAm = float(1.)  # A size parameter in meters
        self.telSizeBm = float(1.)  # B size parameter in meters
        self.etaA = 1. # antenna efficiency (range 0 to 1)
        self.etaB = 1. # efficiency main beam (range 0 to 1)
        self.refChan = refChan
        self.version = str("1.0.1")
        self.polA = str("X")        # polariation of A ydata: X, Y, R, L,
        self.polB = str("Y")        # polariation of B ydata: X, Y, R, L,
        self.polAngle = float(0.0)  # orientation of polariation of A
        self.frame = str("TOPO")    # reference frame (LSR, BARY, TOPO)
# compute coordinates from az,el location and date+time all angles in degrees
        self.ra = float(0.0)        # degrees, ie 12h => 180deg
        self.dec = float(0.0)
        self.gallon = float(0.0)
        self.gallat = float(0.0)
        self.az_sun = float(0.0)
        self.altsun = float(0.0)
        self.epoch = str("2000")
# finally the data
        self.xdata = xdata
        self.ydataA = ydataA
        self.ydataB = ydataB
        self.nChan = len(ydataA)

    def __str__(self):
        """
        Define a spectrum summary string
        """
        secs = self.durationSec
        return "({0}, {1}, {2})".format(self.site, self.utc, str(secs))

    def azel2radec(self):
        """
        Compute the ra,dec (J2000) from Az,El location and time
        """
        location = ephem.Observer()
        location.lon = str(self.tellon)
        location.lat = str(self.tellat)
        location.elevation = self.telelev
        strnow = self.utc.isoformat()
        # convert Time string format into value for Observer
        dates = strnow.split('T')
        datestr = dates[0] + ' ' + dates[1]
        location.date = datestr
        # compute Local Sidereal Time
        lst = location.sidereal_time()
        aparts = angles.phmsdms(str(lst))
        self.lst = angles.sexa2deci(aparts['sign'], *aparts['vals'], todeg=True)
        ## Must set the date before calculating ra, dec!!!
        # compute apparent RA,DEC for date of observations
        ra_a, dec_a = location.radec_of(str(self.telaz), str(self.telel))
        fmt = 'Date   = %s,  LST = %s, %f (%f, %f)'
#        print fmt % (datestr, lst, self.lst, self.telaz, self.telel)
        radec = ephem.Equatorial(ra_a, dec_a, epoch=datestr)
#        print 'Ra,Dec %s,%s for %s' % (radec.ra, radec.dec, radec.epoch)
        radec2000 = ephem.Equatorial( radec, epoch=ephem.J2000)
#        print 'Ra,Dec %s,%s for %s' % (radec2000.ra, radec2000.dec, radec2000.epoch)
        # Hours
        aparts = angles.phmsdms(str(radec2000.ra))
        self.ra = angles.sexa2deci(aparts['sign'], *aparts['vals'], todeg=True)
        # to convert to dec degrees need to replace on : with d
        aparts = angles.phmsdms(str(radec2000.dec))
        self.dec = angles.sexa2deci(aparts['sign'], *aparts['vals'])
        self.epoch = "2000"
        gal = ephem.Galactic(radec2000)
        aparts = angles.phmsdms(str(gal.lon))
        self.gallon = angles.sexa2deci(aparts['sign'], *aparts['vals'])
        aparts = angles.phmsdms(str(gal.lat))
        self.gallat = angles.sexa2deci(aparts['sign'], *aparts['vals'])
        sun = ephem.Sun(location)
        aparts = angles.phmsdms(str(sun.az))
        self.az_sun = angles.sexa2deci(aparts['sign'], *aparts['vals'])
        aparts = angles.phmsdms(str(sun.alt))
        self.altsun = angles.sexa2deci(aparts['sign'], *aparts['vals'])
#        print 'sun az,el: %s,%s -> %f,%f' % (sun.az, sun.alt, self.az_sun, self.altsun)

##################################################
#
    def write_ascii_file(self, dirname, outname):
        """
        Write ascii file containing astronomy data
        """
    # need the current time to update coordiantes
        now = self.utc
        print "File %4d: %s (%d)" % (self.writecount, outname, self.count)
        fullname = dirname + outname
        outfile = open(fullname, 'w')
        outfile.write('# File: ' + outname + '\n')
        if len(str(self.noteA)) > 1:
            outfile.write('# ' + str(self.noteA) + '\n')
        if len(str(self.noteB)) > 1:
            outfile.write('# ' + str(self.noteB) + '\n')
        gainstr = ''
        ngains = len(self.gains)
        for iii in range(ngains-1):
            gainstr = gainstr + str(self.gains[iii]) + '; '
        gainstr = gainstr + str(self.gains[ngains-1])
        outline = '# GAINS     = ' + gainstr + '\n'
        outfile.write(outline)
        outline = '# Count     = ' + str(self.count) + '\n'
        outfile.write(outline)
        outline = '# CenterFreq= ' + str(self.centerFreqHz) + '\n'
        outfile.write(outline)
        outline = '# Bandwidth = '  + str(self.bandwidthHz) + '\n'
        outfile.write(outline)
        outline = '# Duration  = '  + str(self.durationSec) + '\n'
        outfile.write(outline)
        outline = '# DeltaX    = '  + str(self.deltaFreq) + '\n'
        outfile.write(outline)
        nChan = len(self.ydataA)
        outline = '# NCHAN     = '  + str(nChan) + '\n'
        outfile.write(outline)
        strnow = now.isoformat()
        dates = strnow.split('T')
        datestr = dates[0] + ' ' + dates[1]
        outline = '# UTC       = '  + datestr + '\n'
        outfile.write(outline)
        lststr = angles.fmt_angle(self.lst/15., s1=":", s2=":", pre=3)  # convert to hours
        outline = '# LST       = '  + lststr[1:] + '\n'
        outfile.write(outline)
        outline = '# AZ        = '  + str(self.telaz) + '\n'
        outfile.write(outline)
        outline = '# EL        = '  + str(self.telel) + '\n'
        outfile.write(outline)
        anglestr = angles.fmt_angle(float(self.tellon), s1=":", s2=":")
        outline = '# TELLON    = '  + anglestr + '\n'
        outfile.write(outline)
        anglestr = angles.fmt_angle(float(self.tellat), s1=":", s2=":")
        outline = '# TELLAT    = '  + anglestr + '\n'
        outfile.write(outline)
        rastr = angles.fmt_angle(self.ra/15., s1=":", s2=":", pre=3) # convert to hours
        outline = '# RA        = '  + rastr[1:] + '\n'
        outfile.write(outline)
        decstr = angles.fmt_angle(self.dec, s1=":", s2=":")
        outline = '# DEC       = '  + decstr + '\n'
        outfile.write(outline)
        lonstr = angles.fmt_angle(self.gallon, s1=":", s2=":", pre=2)
        outline = '# GALLON    = '  + lonstr[1:] + '\n'
        outfile.write(outline)
        latstr = angles.fmt_angle(self.gallat, s1=":", s2=":", pre=2)
        outline = '# GALLAT    = '  + latstr + '\n'
        outfile.write(outline)
        altstr = angles.fmt_angle(self.altsun, s1=":", s2=":", pre=1)
        outline = '# ALT_SUN   = '  + altstr + '\n'
        outfile.write(outline)
        az_str = angles.fmt_angle(self.az_sun, s1=":", s2=":", pre=1)
        outline = '# AZ_SUN    = '  + az_str + '\n'
        outfile.write(outline)
        outline = '# AST_VERS  = '  + str("02.01") + '\n'
        outfile.write(outline)

        dx = self.bandwidthHz/float(self.nChan)
        x = self.centerFreqHz - (self.bandwidthHz/2.) + (dx/2.)
        yv = self.ydataA
        leny = len(yv)
        for i in range(min(self.nChan,leny)):
            outline = str(i).zfill(4) + ' ' + str(long(x)) + ' ' + str(yv[i]) + '\n'
            outfile.write(outline)
            x = x + dx
        del outline
        outfile.close()

    def write_ascii_ast(self, dirname):
        """
        Write ascii file containing astronomy data
        File name is based on time of observation
        """
        now = self.utc
        strnow = now.isoformat()
        datestr = strnow.split('.')
        # distinguish hot load and regular observations
        if self.telel > 0:
            outname = datestr[0] + '.ast'
        else:
            outname = datestr[0] + '.hot'
        outname = outname.replace(":", "_")
        self.write_ascii_file(dirname, outname)

    def read_spec_ast(self, fullname):
        """
        Read an ascii radio Spectrum file and return a radioSpectrum object
        """

        # Read the file.
        f2 = open(fullname, 'r')
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
            linecount = linecount + 1
# if a very short or blank line
            if len(line) < 3:
                continue
            if linecount == 2:
                self.noteA = line[2:].replace('\n', '')
# if a comment or parameter line, decode value
            if line[0] == '#':
# parse keywords as upper case: ie Ra == RA
                parts[1] = parts[1].upper()
                if parts[1] == 'UTC':
                    timefmt = "%Y-%m-%d %H:%M:%S.%f"
                    utc = datetime.datetime.strptime(parts[3] + " " + parts[4], timefmt)
                    self.utc = utc
                if parts[1] == 'CENTERFREQ':
                    self.centerFreqHz = float(parts[3])
                if parts[1] == 'CENTERFREQ=':
                    self.centerFreqHz = float(parts[2])
                if parts[1] == 'BANDWIDTH':
                    self.bandwidthHz = float(parts[3])
                if parts[1] == 'DELTAX':
                    self.deltaFreq = float(parts[3])
                if parts[1] == 'LST':
                    lstparts = angles.phmsdms(parts[3])
                    x = angles.sexa2deci(lstparts['sign'], *lstparts['vals'])
                    self.lst = x*15. # convert back to degrees
#                    print parts[3], x
                if parts[1] == 'AZ':
                    self.telaz = degree2float(parts[3], parts[1])
                if parts[1] == 'EL':
                    self.telel = degree2float(parts[3], parts[1])
                if parts[1] == 'COUNT':
                    self.count = int(parts[3])
                if parts[1] == 'NCHAN':
                    self.nChan = int(parts[3])
                if parts[1] == 'REFCHAN':
                    self.refChan = float(parts[3])
                if parts[1] == 'ETAA':
                    self.etaA = float(parts[3])
                if parts[1] == 'ETAB':
                    self.etaB = float(parts[3])
                if parts[1] == 'POLANGLE':
                    self.polAngle = float(parts[3])
                if parts[1] == 'LNA' or parts[1] == 'GAINS':  # get one or more gains separated by ';'
                    gains = []
                    for jjj in range(3, len(parts)):
                        gainstr = parts[jjj].replace(';', ' ')
                        gainstr = gainstr.replace(',', ' ')
                        moreparts = gainstr.split()
                        for kkk in range(len(moreparts)):
                            gains.append(float(moreparts[kkk]))
                    self.gains = numpy.array(gains)
                if parts[1] == 'LNA=' or parts[1] == 'GAINS=':  # get one or more gains separated by ';'
                    gains = []
                    for jjj in range(2, len(parts)):
                        gainstr = parts[jjj].replace(';', ' ')
                        gainstr = gainstr.replace(',', ' ')
                        moreparts = gainstr.split()
                        for kkk in range(len(moreparts)):
                            gains.append(float(moreparts[kkk]))
                    self.gains = numpy.array(gains)
                if parts[1] == 'DURATION':
                    self.durationSec = time2float(parts[3], parts[1])
                if parts[1] == 'LON' or parts[1] == 'GALLON':
                    aparts = angles.phmsdms(parts[3])
                    x = angles.sexa2deci(aparts['sign'], *aparts['vals'])
                    self.gallon = x
                if parts[1] == 'LAT' or parts[1] == 'GALLAT':
                    aparts = angles.phmsdms(parts[3])
                    x = angles.sexa2deci(aparts['sign'], *aparts['vals'])
                    self.gallat = x
# if parse telescope geographic latitude and longitude into float
                if parts[1] == 'TELLON':
                    self.tellon = degree2float(parts[3], parts[1])
                if parts[1] == 'TELLAT':
                    self.tellat = degree2float(parts[3], parts[1])
# parse ra, dec into float
                if parts[1] == 'RA':
                    aparts = angles.phmsdms(parts[3])
                    x = angles.sexa2deci(aparts['sign'], *aparts['vals'])
                    self.ra = x * 15. # convert back to degrees
#                    print 'RA', parts[3], aparts, x
                if parts[1] == 'DEC':
                    aparts = angles.phmsdms(parts[3])
                    x = angles.sexa2deci(aparts['sign'], *aparts['vals'])
                    self.dec = x
# if sun coordinates into a float
                if parts[1] == 'ALT_SUN':
                    aparts = angles.phmsdms(parts[3])
                    x = angles.sexa2deci(aparts['sign'], *aparts['vals'])
                    self.altsun = x
                if parts[1] == 'AZ_SUN':
                    aparts = angles.phmsdms(parts[3])
                    x = angles.sexa2deci(aparts['sign'], *aparts['vals'])
                    self.az_sun = x
                continue
# sometimes there are user comments in the top few lines; ignore
            if linecount < 5:
                continue
# start data processing
            datacount = datacount+1
            p = line.split()
            x1.append(float(p[1]))
            y1.append(float(p[2]))

# at this point all data and header keywords are read
#        self.xdata = numpy.array(x1)*1.e6  # convert to Hz
        self.xdata = numpy.array(x1)  # convert to Hz
        self.ydataA = numpy.array(y1) #scale to a reasonable value for plotting
        ndata = len(self.xdata)
        if self.nChan != ndata:
            print "File header Miss-match and number of channels in data"
            print ": %f != %f" % (self.nChan, ndata)
            self.nChan = int(ndata)
        return

# HISTORY
# 16AUG17 GIL separate out modules and ephemeris
# 15JUN04 GIL recording working, but having trouble with display
# 15JUN02 GIL return coordinates so that the data can be understood
# 15MAY29 GIL get averaging functioning and also periodic updates
# 15MAY16 GIL start creating averages
# 15MAY15 GIL try to reduce the cpu load
# 15MAY05 GIL have recording working now.  still too much CPU usage.
# 15MAY04 GIL found that the program was taking too much CPU
#             performance was fine once sleep of 10ms was added
# 15MAY01 GIL Initial version


#Example file is below;  Comments (#) are part of the file, except for the data columns
# File: 2016-05-06T23:13:39.ast
# Time range: 2016-05-06T23:13:39.242442 - 2016-05-06T23:11:39.231024
# P162LN+P103LN+3 LNA4ALL +P103LN+ testhorn + 1445filter
# LNA= 12.5;
# Count     = 1798
# CenterFreq= 1420550000.0
# Bandwidth = 2500000.0
# DeltaX    = 2441.40625
# NCHAN     = 1024
# UTC       = 2016-05-06 23:13:39.242442
# LST       = 21:16:41.03
# AZ        = 230.0
# EL        = -90.0
# RA        = 9:16:03.40
# DEC       = -41:58:46.3
# LON       = 265:45:28.6
# LAT       = 4:51:32.4
# ALT_SUN   = 14:45:23.7
# AZ_SUN    = 80:22:20.0
# .... Data starts here .... (NO # in actual file
###"""
###0000 1419300000 7.89901056399e-08
###0001 1419302441 7.70544758831e-08
###"""
