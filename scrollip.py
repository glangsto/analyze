#!/usr/bin/python
# python to get my ip and send an email to megajansky
#HISTORY
# 21Jan16 GIL exit if can not reach scroll hat
# 20Dec23 GIL reduce text scolled
# 20Dec19 GIL add stop option
# 20Dec15 GIL initial version
# 20Nov10 GIL try to readh web ip
# 20Nov03 GIL initial version

import os
import sys
import time
try:
    from sense_hat import SenseHat
except:
    sys.exit()

# configure sense hat
try:
    sense = SenseHat()
except:
    print("Can not reach Scroll Hat")
    sys.exit()
try:
    sense.set_rotation(0)
except:
    print("Can not set hat rotation")
    sys.exit()
    
red = (155, 0, 0)
green = (0, 55, 0)
blue = (50, 50, 255)
purple = (155, 0, 155)
sense.show_message("SCROLLIP -", text_colour=green)
sense.show_message(" Langston", text_colour=green)
sense.show_message(" NSF", text_colour=blue)

waitcount = 0

# communicate with this code that runs as root, via /tmp
# directory
mynewevent = "/tmp/mynewevent"
mystop = "/tmp/stop"

while True: 

    # when the temp file exists, then flash an event
    if os.path.exists(mynewevent):
        os.system("sudo rm -f %s" % (mynewevent))
        sense.show_message("*", text_colour=red)
#        time.sleep(1)
    else:
        time.sleep(1)

    if os.path.exists(mystop):
        os.system("sudo rm -f %s" % (mystop))
        sys.exit()

# only print IP every so often, but watch for events        
    if waitcount < 1:
        myiptemp = "/tmp/myip"
        os.system( "/home/pi/bin/myip | head -1 > %s" % (myiptemp))
        # read all the ip addresses
        with open(myiptemp, "r") as f:
            for line in f:
                ipparts = str(line)

        f.close()
        os.system( "sudo rm -f %s" % (myiptemp))

        parts = ipparts.split('\n')
        myip = parts[0]
        mytext = "IP: " + myip
        sense.show_message(mytext, text_colour=purple)

        
    waitcount = waitcount + 1

    # now time to flash the ip
    if waitcount > 10:
        waitcount = 0
        # temperature
        temp = sense.get_temperature()
        tempstr = "T: %5.1f C" % (float(temp))
        sense.show_message(tempstr, text_colour=red)

        # humidity
        humidity = sense.get_humidity()
        humidstr = "H: %4.0f %%" % (float(humidity))
        sense.show_message(humidstr, text_colour=blue)
        
        # humidity
        pressure = sense.get_pressure()
        pressstr = "P: %5.0f mbar" % (float(pressure))
        sense.show_message(pressstr, text_colour=purple)
        

        



