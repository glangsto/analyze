#!/usr/bin/env python
"""
Create a full path to a new directory, creating parent directories if needed
Usage: python fullpath.py my/new/directory
"""
import os
import sys
from pathlib import Path

def main():
    """ Parse arguments, then fix whitespace in the given file """
    if len(sys.argv) == 2:
        pathname = sys.argv[1]
    else:
        print("Invalid arguments. Usage: python fullpath.py my/new/direcotry")
        sys.exit(1)
    fullpath(pathname)


def fullpath(pathname):
    """
    Create all directories leading up to the new directory
    """
    newpath = Path(pathname)

    newpath.mkdir( parents=True, exist_ok=True)

    print(f"Path created: {newpath}")

if __name__ == "__main__":
    main()
