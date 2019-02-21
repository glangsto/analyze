# Script to convert a Cologne Molecular Spectroscopy Database
# Search output file into a input file for the rt_idl tmcLineSum
# 
# HISTORY
# 17OCT12 GIL add sum of all line strenghts used 
# 14FEB26 GIL weight by 1/line-strength; keep transition numbers < 100
# 14JAN07 GIL update for AAS meeting
# 13NOV12 GIL update weights estimate to SNR^2
# 13NOV10 GIL add in weighting estimate list
# 13NOV08 GIL first completed version, used to make many lists

import sys
import math

prefix = 'c-C3HCN'
prefix = 'c6h'
#
# get the command line arguments
args = sys.argv
nargs = len(args)
if nargs <= 1:
    print "cdms takes an Cologne Molecular Spectrocopy Database"
    print "Model Line list and converts it into a set of IDL variables"
    print "These line frequency, and intensity lists are used to search"
    print "Spectral serveys for lines"
    print ""
    print "Usage:"
    print "cdms [-r MinFreqMHz MaxFreqMHz] <lineList> [<fraction>]"
    print "where"
    print "<lineList>.txt is the input file"
    print "               and   <lineList>.pro is the output file"
    print "<fraction>  fractional minimum of strongest line strength to keep"
    print "<fraction>  range 0.0 to 1.0, default 0.0"
    print ""
    print "Glen Langston (NSF) 2013 December 10"
    exit
# if first argument is not flagging the frequency range

minFreq = 0.
maxFreq = 300000.

# the first argument is the molecule input file list
prefix = args[1]
if prefix.lower() == '-r':
    minFreq = float(args[2])
    maxFreq = float(args[3])
    prefix = args[4]
    print "Limiting input lines to between "+str(minFreq)+" and "+str(maxFreq)+" MHz"
print 'Prefix: ', prefix
molecule=prefix.upper()
fraction = 0.0   # fraction limit 0 == all lines
if (nargs == 3) or (nargs == 6):
    fraction = float(args[nargs-1])

if (fraction > .99) or (fraction < 0.):
    fraction = 0.
#
print 'Fraction: ', fraction

infile = prefix + '.txt'
outfile = prefix + '.pro'

f   = open( infile, 'r')
out = open( outfile, 'w')
print >> out, "plotTitle = '%s' " % (molecule)
print >> out, "plotFile = '%s' " % (prefix)
print >> out, "inFreqs = [  $\\ "

nIn=0                # get min, max intensity, count  lines
nLabel = 0
nFreq = 0
nWeight = 0
K = 0.0
maxK = 0.0
sumK = 0.0
minK = 10000.0
# for all known molecular lines in file 
for line in f:
    # print line,
    words = line.split()
    nwords = len(words)
    if (nwords > 2):
        if (words[0] == "!"):
            print 'Comments: ',line
            continue
        nIn = nIn + 1
        nu = float( words[0])
        freq = float(nu)
        if (freq < minFreq) or (freq > maxFreq):
            continue
        intens = words[2]
        K = math.pow( 10., float(intens))
        if K > maxK:
            maxK = K
        if K < minK:
            minK = K
        sumK = sumK + K
f.close()

# now min and max model intensities known
f = open( infile, 'r')
nout = 0
for line in f:
    # print line,
    words = line.split()
    nwords = len(words)
    if nwords > 2:
        if (words[0] == "!"):
            continue
        nu = words[0]
        freq = float( nu)
        if (freq < minFreq) or (freq > maxFreq):
            continue
        name = words[nwords-1]
        dnu = words[1]
        intens = words[2]
        K = math.pow( 10., float(intens))
        snr2 = pow(K/maxK,1.) # weighting squared ratio
        if snr2 > fraction:
            buf = " , %s $\\ ; +/-%s %s %s" % (nu, dnu, name, intens)
            if (nout == 0):
                buf = "   %s $\\ ; +/-%s %s %s" % (nu, dnu, name, intens)
            print >> out, buf
            nout = nout + 1
            nFreq = nFreq + 1
print >> out, " ] "
print >> out, " "

f.close()
print "Line Intensity Range: ", minK, " to ", maxK, " (", nout, " )"
print "Sum of all Line Strengths: ", sumK

f = open( infile, 'r')

noLabel = 1

nout= 0
print >> out, "inLabels = [ $\\ "

for line in f:
    # print line,
    words = line.split()
    nwords = len(words)
    if nwords > 2:
        if (words[0] == "!"):
            continue
        nu = words[0]
        freq = float(nu)
        if (freq < minFreq) or (freq > maxFreq):
            continue
        intens = words[2]
        if (noLabel < 1):
            trans1 = int(words[nwords-4])
            trans2 = int(words[nwords-3])
            trans3 = int(words[nwords-2])
            name = words[nwords-1]
        else:
            trans1 = int(words[nwords-4])
            trans2 = int(words[nwords-3])
            trans3 = int(words[nwords-2])
            name = ""
        trans1 = trans1 % 100
        trans102 = trans2 % 100
        trans103 = trans3 % 100
        # if transtions > 100, its probably two transitions
        if (trans3 > trans103):
            trans1 = trans2 
            trans2 = int( trans3/100)
            trans3 = trans103
            trans102 = trans2 % 100
        if (trans2 > trans102):
            trans1 = trans2 / 100
            trans2 = trans102
        trans1 = trans1 % 100
        K = math.pow( 10., float(intens))
        snr2 = pow(K/maxK,1.) # weighting squared ratio
        if snr2 > fraction:
            buf = "  '%s %s %s %s'  $\\ ; %s " % (name, trans1, trans2, trans3, nu)
            if (nout > 0): 
                buf = ", '%s %s %s %s'  $\\ ; %s " % (name, trans1, trans2, trans3, nu)
            print >> out, buf
            nout = nout + 1
            nLabel = nLabel + 1
print >> out, " ] "
print >> out, " "

f.close()
f   = open( infile, 'r')
print >> out, "weights = [ $\\ "

nout = 0

for line in f:
    # print line,
    words = line.split()
    nwords = len(words)
    if nwords > 2:
        if (words[0] == "!"):
            continue
        nu = words[0]
        freq = float( nu)
        if (freq < minFreq) or (freq > maxFreq):
            continue
        intens = words[2]
        trans1 = words[nwords-4]
        trans2 = words[nwords-3]
        trans3 = words[nwords-2]
        name = words[nwords-1]
        K = math.pow( 10., float(intens))
#        snr2 = pow(K/maxK,1.) # weighting squared ratio
# temporary test
        snr2 = K/maxK
        if snr2 > fraction:
            buf = "  %8.5f $\\ ; %s %s %s %s " % (snr2, trans1, trans2, trans3, nu)
            if (nout > 0): 
                buf = ", %8.5f $\\ ; %s %s %s %s " % (snr2, trans1, trans2, trans3, nu)
            print >> out, buf
            nout = nout + 1
            nWeight = nWeight + 1
print >> out, " ] "
print >> out, " "

f.close()
print >> out, "Nout = ",nout
idlprint = 'print, "Read frequencies for ",Nout," Lines'
print >> out, idlprint
out.close()

print 'Wrote frequencies for ', nout, ' of ', nIn,' Lines'
print 'NFreq = ', nFreq, '; nLabel = ', nLabel,'; nWeight = ', nWeight
