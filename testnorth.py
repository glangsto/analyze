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
    if waitcount > 2:
        waitcount = 0
        # temperature
        sense.set_imu_config(True, True, True)
        orientation = sense.get_orientation()
        pitch = float(orientation['pitch'])
        roll = float(orientation['roll'])
        yaw = float(orientation['yaw'])
#        print("p: {pitch}, r: {roll}, y: {yaw}".format(**orientation))

        # humidity
        sense.set_imu_config(True, False, False)
        north = float(sense.get_compass())
        print("p: %7.1f, r: %7.1f, y: %7.1f, North: %7.1f" % \
              (pitch, roll, yaw, north))
    else:
        time.sleep(1)

    waitcount = waitcount + 1

        



