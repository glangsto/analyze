# Save as server.py 
# Message Receiver
import os
from socket import *
import datetime

host = ""
port = 13000
buf = 1024
addr = (host, port)
UDPSock = socket(AF_INET, SOCK_DGRAM)
UDPSock.bind(addr)
print "Waiting to receive messages..."
currentDT = datetime.datetime.now()
print "Starting Server now: ", str(currentDT)

while True:
    (data, addr) = UDPSock.recvfrom(buf)
    currentDT = datetime.datetime.now()
    print "Received message: " + data
    print "Received at time:", str(currentDT)
    if data == "exit":
        break
UDPSock.close()
os._exit(0)
