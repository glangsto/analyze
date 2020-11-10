# python script to get the location of the IP server and
# print out the latitude and longitude in a clearer format
# HISTORY
# 20Nov03 GIL Clean up dependances.

import re
import json
import os
from urllib2 import urlopen

try:
    data = str(urlopen('http://checkip.dyndns.com/',timeout=10).read())

    IP = re.compile(r'(\d+.\d+.\d+.\d+)').search(data).group(1)
    url = 'http://ipinfo.io/' + IP + '/json'
    response = urlopen(url)
    data = json.load(response)

    org=data['org']
    city = data['city']
    country=data['country']
    region=data['region']
    loc = data['loc']
    locs = loc.split(',')
    flat = float( locs[0])
    flon = float( locs[1])

    print 'Your IP details\n '
    print 'IP : {4} \nRegion : {1} \nCountry : {2} \nCity : {3} \nOrg : {0}'.format(org,region,country,city,IP)

except:
    flat = 0.
    flon = 0.
    IP = "127.0.0.1"
    
webiptemp = "/tmp/webip"
os.system( "sudo rm %s" % (webiptemp))
os.system( "echo %s > %s" % (IP, webiptemp))

# keep latitude and longitude in float format
lat = flat
lon = flon

# convert to degrees minutes seconds
# get sign of latitude
if lat >= 0:
    pmlat = '+'
else:
    pmlat = '-'
    lat = - lat

# get sign of longitude
if lon >= 0:
    pmlon = '+'
else:
    pmlon = '-'
    lon = - lon

# get degrees part
dlat = int(lat)
dlon = int(lon)

# start on minutes part
lat = (lat - dlat)*60.
lon = (lon - dlon)*60.

mlat = int(lat)
mlon = int(lon)

# start on seconds part
slat = (lat - mlat)*60.
slon = (lon - mlon)*60.

# show results
print('Latitude : %9.5f' % (flat))
print('Latitude : %s%02d:%02d:%05.2f' % (pmlat,dlat,mlat,slat))
print('Longitude: %9.5f' % (flon))
print('Longitude: %s%02d:%02d:%05.2f' % (pmlon,dlon,mlon,slon))


