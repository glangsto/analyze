#!/usr/bin/env python3
"""
strip trailing whitespace from file
usage: stripspace.py <file>
"""

import sys

if len(sys.argv[1:]) != 1:
  sys.exit(__doc__)

content = ''
outsize = 0
inp = outp = sys.argv[1]
with open(inp, 'rb') as infile:
  content = infile.read()
with open(outp, 'wb') as output:
  for line in content.splitlines():
    newline = line.rstrip()
    outsize += len(newline) + 1
    output.write(newline + '\n')

print("Done. Stripped %s bytes." % (len(content)-outsize))
