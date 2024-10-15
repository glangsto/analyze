#!/usr/bin/env python
"""
Fix trailing whitespace and line endings (to Unix) in a file.
Usage: python fix_whitespace.py foo.py
"""

import os
import sys


def main():
    """ Parse arguments, then fix whitespace in the given file """
    if len(sys.argv) == 2:
        fname = sys.argv[1]
        if not os.path.exists(fname):
            print("Python file not found: %s" % sys.argv[1])
            sys.exit(1)
    else:
        print("Invalid arguments. Usage: python fix_whitespace.py foo.py")
        sys.exit(1)
    fix_whitespace(fname)


def fix_whitespace(fname):
    """ Fix whitespace in a file """
    with open(fname, "rb") as fo:
        original_contents = fo.read()
    # "r"
    with open(fname, "r") as fo:
        contents = fo.read()
    lines = contents.split("\n")
    fixed = 0
    nlbytes = str.encode("\n")
    for k, line in enumerate(lines):
        new_line = line.rstrip()
        if len(line) != len(new_line):
            lines[k] = new_line
            fixed += 1
    with open(fname, "wb") as fo:
##        for aline in lines:
#            alinebytes = str.encode(aline)
        fo.write(nlbytes.join(line))
    if fixed or contents != original_contents:
        print("************* %s" % os.path.basename(fname))
    if fixed:
        slines = "lines" if fixed > 1 else "line"
        print("Fixed trailing whitespace on %d %s" \
              % (fixed, slines))
    if contents != original_contents:
        print("Fixed line endings to Unix (\\n)")


if __name__ == "__main__":
    main()
