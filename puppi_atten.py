"""
Routines for setting (manually or automatically) the PUPPI attenuators
"""
import numpy as np
import socket
import guppi.script as gup
import time
import logging

def autoPuppiAtten(goal_rms=20.0, max_tries = 3):
    tries = 0
    attenA,attenB = getPuppiAttens()
    while tries < max_tries:
        gup.arm()
        alev = gup.get_adc_samples(fpga=3).std()
        blev = gup.get_adc_samples(fpga=1).std()
        
        corrA = 20*np.log10(alev/goal_rms)
        corrB = 20*np.log10(blev/goal_rms)
        logging.info( "Currently at: A atten=%.1f, rms=%.1f, error=%.1fdB  |  B atten=%.1f, rms=%.1f, error=%.1fdB" % (attenA, alev, corrA, attenB, blev, corrB))
        if np.abs(corrA) < 2.0 and np.abs(corrB) < 2.0:
            print "Both channels within 2 dB of goal; finished"
            logging.info("auto atten successful with final values %.1f %.1f" % (attenA,attenB))
            break
        if np.abs(corrA) > 2.0:
            newA = atten_range(attenA + corrA)
        else:
            newA = attenA
        if np.abs(corrB) > 2.0:
            newB = atten_range(attenB + corrB)
        else:
            newB = attenB
        setPuppiAttens(newA, newB)
        attenA = newA
        attenB = newB
        tries += 1
    else:
        print "Failed to auto level"
        

def getPuppiAttens(ip_addr = 'ifamp', ip_port = 6001):
    s = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
    s.settimeout( int(5) )
    try:
        s.connect((ip_addr, int(ip_port))) 
        logging.debug("-- socket: opening socket %s/%d " %(ip_addr, ip_port))
    except:                          # Time out
        logging.debug("-- socket: no connection to %s/%d " %(ip_addr, ip_port))
        sys.exit()
    resp = send_recv(s,'ATN?')
    if resp.find('atnm') != 0:
        print "error parsing attenuator response"
    a = float(resp[4:6])/2.0
    b = float(resp[6:8])/2.0
    s.close()
    return a,b

def setPuppiAttens(a,b, ip_addr = 'ifamp', ip_port = 6001):
    s = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
    s.settimeout( int(5) )
    try:
        s.connect((ip_addr, int(ip_port))) 
        logging.debug("-- socket: opening socket %s/%d " %(ip_addr, ip_port))
    except:                          # Time out
        logging.debug("-- socket: no connection to %s/%d " %(ip_addr, ip_port))
        sys.exit()
        
    udc_command = "ATNA%2.2d" % (a * 2.0)
    resp = send_recv(s,udc_command)
    if (resp != "atnok\r"):
        print "Error setting atten_a: %s" % (resp)
        
    udc_command = "ATNB%2.2d" % (b * 2.0)
    resp = send_recv(s,udc_command)
    if (resp != "atnok\r"):
        print "Error setting atten_a: %s" % (resp)
    
    s.close()

    
def atten_range(db):
    if (db>15.5):
        db=15.5
    if (db<0.0):
        db=0.0
    return db

def send_recv(sock,message):
    BUFFER_SIZE = 1024
    try:
        sock.send(message + '\r')
    except:
        print "-- error sending message '%s'" % (message)
    time.sleep(0.5)
    try:
        response = sock.recv(BUFFER_SIZE)
        return response
    except:
        print "-- timeout receiving response to '%s'" % (message)
        return None
