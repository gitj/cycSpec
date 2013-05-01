import socket
import time
import numpy as np

def getUdpSocket(port=61000):
    s = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    s.bind(('',port))
    return s
    
def rxData(s):
    lasth = 0
    while True:
        d = s.recv(10000)
        if len(d) != 8200:
            print time.ctime(), (" %.3f runt: %016x %d" % (np.fmod(time.time(),1),lasth,len(d)))
        h = np.fromstring(d[-8:],dtype='>u8')[0]
        if h - lasth != 8192:
            print time.ctime(), (" %.3f skip: %016x %016x %d" % (np.fmod(time.time(),1),h,lasth,h-lasth))
        lasth = h