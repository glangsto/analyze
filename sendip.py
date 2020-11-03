# python to get my ip and send an email to megajansky
# HISTORY
# 20Nov03 GIL initial version

import os

myiptemp = "/tmp/myip"
os.system( '/home/pi/bin/myip | head -1 > /tmp/myip')

# prepare for failure to read parts of ip
ipparts = [ "1", "2"]

# read all the ip addresses
with open(myiptemp, "r") as f:
    for line in f:
        ipparts = str(line)

f.close()

parts = ipparts.split('\n')
myip = parts[0]
#print ("My IP: %s" % (myip))
# 
ipparts = "/tmp/myip." + myip

parts = ipparts.split('\n')
ipfile = parts[0]
#now write ip to a temp file
mailtemp = "/tmp/myemail"
os.system( "echo 'FROM: root' > %s" % (mailtemp))
os.system( "echo 'TO: root' >> %s" % (mailtemp))
os.system( "echo 'Subject: Pi Booting: %s' >> %s" % (myip,mailtemp))
os.system( "echo 'Log of Pis booting up: %s' >> %s" % (myip, mailtemp))

# read all the user-suppied horn info into the mail log
os.system( "cat /boot/horn.txt >> %s" % (mailtemp))
# seem to need to wait for network to get working, before getting ip
os.system( "/bin/sleep 15")
# now get the location of the ip address
os.system( "/home/pi/Research/analyze/iplatlon | cat >> %s" % (mailtemp))

os.system( "cat %s | /usr/bin/sendmail -t" % (mailtemp))

