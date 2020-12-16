#!/usr/bin/python
# python to get my ip and send an email to megajansky
#HISTORY
# 20Dec15 GIL initial version
# 20Nov10 GIL try to readh web ip
# 20Nov03 GIL initial version

import os
import sys
import time
from sense_hat import SenseHat

# configure sense hat
sense = SenseHat()
sense.set_rotation(0)
red = (155, 0, 0)
green = (0, 55, 0)
blue = (50, 50, 255)
purple = (155, 0, 155)
sense.show_message("scrollip -- ", text_colour=green)
sense.show_message(" Glen Langston ", text_colour=green)
sense.show_message(" NSF 20 Dec 15 ", text_colour=blue)

waitcount = 0

while True: 

    # when the temp file exists, then flash an event
    mynewevent = "/tmp/mynewevent"
    if os.path.exists(mynewevent):
        os.system("sudo rm -f %s" % (mynewevent))
        sense.show_message("*", text_colour=red)
#        time.sleep(1)
    else:
        time.sleep(1)

# only print IP every so often, but watch for events        
    if waitcount < 1:
        myiptemp = "/tmp/myip"
        os.system( "/home/pi/bin/myip | head -1 > %s" % (myiptemp))
        # read all the ip addresses
        with open(myiptemp, "r") as f:
            for line in f:
                ipparts = str(line)

        f.close()

        parts = ipparts.split('\n')
        myip = parts[0]
        mytext = "IP: " + myip
        sense.show_message(mytext, text_colour=purple)

    waitcount = waitcount + 1

    # now time to flash the ip
    if waitcount > 10:
        waitcount = 0
