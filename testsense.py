#!/usr/bin/python
# python to get my ip and send an email to megajansky
#HISTORY
# 20Dec19 GIL add stop option
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

# communicate with this code that runs as root, via /tmp
# directory
mynewevent = "/tmp/mynewevent"
mystop = "/tmp/stop"
waitcount = 100

while True: 

    # now time to flash the ip
    if waitcount > 10:
        waitcount = 0
        # temperature
        temp = sense.get_temperature()
        tempstr = "Temperature: %5.1f C" % (float(temp))
#        sense.show_message(tempstr, text_color=red)
        print(tempstr)

        # humidity
        humidity = sense.get_humidity()
        humidstr = "Humiidty:    %5.1f %%" % (float(humidity))
        print(humidstr)
#        sense.show_message(humidstr, text_color=blue)
        
        # humidity
        pressure = sense.get_pressure()
        pressstr = "Pressure:    %5.1f Millibars" % (float(pressure))
#        sense.show_message(pressstr, text_color=purple)
        print(pressstr)
    else:
        time.sleep(1)

    waitcount = waitcount + 1

        



