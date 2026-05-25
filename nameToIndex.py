#!/usr/bin/env python
#HISTORY
#26May24 GIL initial version
"""
Find the telescope/cpu index in the file name
Currently only works with names like: pi11-data-26May24
"""
def getIndex(pathname):
    """
    determine if a "pi" is in the name, then find the index
    if no telescopeindex return 0
    """

    parts = pathname.split("pi")
    if len(parts) < 2:
        print("No pi in file name")
        return (0)

    # now should have a name like
    "11-data-26May24"
    name = parts[1]
    nameparts = name.split("-")
    try:
        telescopeIndex = int(nameparts[0])
    except:
        telescopeIndex = 0
        
    return telescopeIndex
                        
