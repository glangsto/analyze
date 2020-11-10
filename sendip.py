# python to get my ip and send an email to megajansky
# HISTORY
# 20Nov10 GIL try to readh web ip
# 20Nov03 GIL initial version

import os

myiptemp = "/tmp/myip"
os.system( "/home/pi/bin/myip | head -1 > %s" % (myiptemp))
    # read all the ip addresses
with open(myiptemp, "r") as f:
    for line in f:
        ipparts = str(line)

f.close()

parts = ipparts.split('\n')
myip = parts[0]

webip = "127.0.0.1"
webiptemp = "/tmp/webip"

# prepare for failure to read parts of ip
ipparts = [ "1", "2"]

try:
    # read all the ip addresses
    with open(webiptemp, "r") as f:
        for line in f:
            ipparts = str(line)

    f.close()

    if len(ipparts) > 0:
        parts = ipparts.split('\n')
        webip = parts[0]
    
except:
    webip = myip
    
print ("My IP: %s, Web IP: %s" % (myip, webip))

#now write ip to a temp file
mailtemp = "/tmp/myemail"
os.system( "echo 'FROM: root' > %s" % (mailtemp))
os.system( "echo 'TO: root' >> %s" % (mailtemp))
os.system( "echo 'Subject: Pi Booting: %s' >> %s" % (webip,mailtemp))
os.system( "echo 'Log of Pis booting up: %s' >> %s" % (webip, mailtemp))

# read all the user-suppied horn info into the mail log
os.system( "cat /boot/horn.txt >> %s" % (mailtemp))
# read all the last IP this pi had
os.system( "cat /boot/lastip.txt >> %s" % (mailtemp))
# seem to need to wait for network to get working, before getting ip
#os.system( "/bin/sleep 15")
# now get the location of the ip address
os.system( "/home/pi/Research/analyze/iplatlon | cat >> %s" % (mailtemp))

os.system( "cat %s | /usr/sbin/sendmail -t" % (mailtemp))

# now clean up for next go
os.system( "sudo rm -f %s" % (myiptemp))
#os.system( "sudo rm -f %s" % (webiptemp))
os.system( "sudo rm -f %s" % (mailtemp))
