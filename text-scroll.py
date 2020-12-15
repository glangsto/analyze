#!/usr/bin/python
# Scroll text on the Pi Sense Hat.
#HISTORY
#20Dec15 GIL initial version
import os
import sys
import time

from sense_hat import SenseHat

nargs = len(sys.argv)
doLoop = False

if nargs < 2:
    mytext = "Hi Pi"
else:
    mytext = sys.argv[nargs-1]
    anArg = sys.argv[1].upper()
    if anArg == "-L":
        doLoop = True


sense = SenseHat()
sense.set_rotation(0)
red = (155, 0, 0)
purple = (155, 0, 155)

while True: 
    sense.show_message(mytext, text_colour=purple)
    time.sleep(1)
    if not doLoop:
        break
