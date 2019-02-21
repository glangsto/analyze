#Python Script to fix the azimuth of an astronomy file
#plot the raw data from the observation
#HISTORY
#18Feb27 GIL initial version

import sys
import statistics
import radioastronomy
import numpy as np
dy = -1.

nargs = len( sys.argv)

if nargs < 2:
    print 'FIXAZ: correct the logged azimuth of an observation'
    print 'FIXAZ <azimuth> filenames'
    print 'where azimuth Azimuth of telescope in degrees.'
    exit()

azimuth = float(sys.argv[1])

names = sys.argv[2:]

n = len(names)

rs = radioastronomy.Spectrum()
for filename in names:

    parts = filename.split('.')
    if parts < 2:
        continue
    nparts = len(parts)
    extension = parts[nparts-1]
    extension = extension.upper()
    # can only process astronomical files
    if (extension != 'AST') and (extension != 'HOT') and (extension != 'CLD'):
        print 'Non-astronomy files can not be read: ',filename,' skipping'
        continue
    rs.read_spec_ast( filename)
    rs.telaz = azimuth
    rs.azel2radec()    # compute ra,dec from az,el

    rs.write_ascii_file( "./", filename)
