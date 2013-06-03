#!/bin/env python
import numpy as np
from matplotlib import pyplot as plt
import os
import socket
import guppi.script as gup
import time
from parfile import psr_par
import sys
from guppi_daq import astro_utils
import observe_utils

useOsf = True # should we do real-time cyclic spectroscopy when obserivng at 327 and 430 MHz?

import logging
logging.basicConfig(filename=('/home/gpu/gjones/logs/nanoObserve.log'),level=logging.DEBUG,
                        format = (':%(levelname)s - %(asctime)s - %(module)s:%(funcName)s:%(lineno)d - %(message)s'))
logger = logging.getLogger()
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(logging.Formatter(':%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(sh)

osfnodes = [1,2,8,9]
puppinodes = [3,5,6,7]



guppiConfigs = {'327' : ('c64', '100'),
                '430' : ('c64', '100'),
                'L' : ('c512', '800'),
                'S' : ('c512', '800'),
                }

rcvrConfigs = {"L": dict(
                    freq = 1380.0,
                    nchan = 512,
                    bw = -800.0,
                    pol = "LIN",
                    rcvr = "L-wide"),
               "S": dict(
                    freq = 2030.0,    # with low-band filter, better aliasing
                    nchan = 512,
                    bw = -800.0,
                    pol = "LIN",
                    rcvr = "S-wide"),
               "327": dict(
                    freq = 327.0,
                    nchan = 64,
                    bw = 100.0,
                    pol = "LIN",
                    rcvr = "327"),
               "430": dict (
                    freq = 430.0,
                    nchan = 64,
                    bw = 100.0,
                    pol = "CIRC",
                    rcvr = "430")
               }


def observe(session,endAt,projid, calSession, useOsf = useOsf,test=False):
    if useOsf:
        logging.info("Checking real time cyclic spectroscopy hardware")
        try:
            import osf
            osf.getServers(osfnodes)
        except:
            logging.exception("Error with setting up cyclic spectroscopy hardware.")
            raw_input("Press enter to proceed without cyclic spectroscopy hardware, or press ctrl-c to stop the script. Then re-run the script after fixing the problem.")
            useOsf = False
            
        try:
            osf.checkRoach()
            osf.stopData()
            osf.killAll(osfnodes)
        except:
            logging.exception("Found ROACH not in expected state; trying to reprogram")
            try:
                osf.setAdcClock(512.0)
                osf.init()
                osf.checkRoach()
            except:
                logging.exception("Failed to reinitialize ROACH!")
                raw_input("Press enter to proceed without cyclic spectroscopy hardware, or press ctrl-c to stop the script. Then re-run the script after fixing the problem.")
                useOsf = False

        if useOsf:
            print "\n\n" 
            logging.info("Cyclic spectroscopy hardware is configured and ready to go!")
            
    lastMode = None # most recent puppi mode, should be None normally

    if calSession:
        lastBand = None
        logging.info("Doing flux calibration session")
        for (source,parfile,band) in calSession:
            logging.info("\n\n*** Slew to source: "+ source)
            if lastBand != band:
                print "\n*** Select band:", band
                if raw_input("Press enter when finished. (Type 'skip' to skip configuring puppi)") != 'skip':
                    mode,clk = guppiConfigs[band]
                    if mode != lastMode:
                        logging.info("Configuring puppi with mode " + mode +" and clock " + clk)
                        gup.mode(mode)
                        gup.reset(clk)
                        gup.arm()
                        lastMode = mode
                    else:
                        logging.info("previous PUPPI mode was the same as current, so not reconfiguring. Mode is: " + mode)
        
                
                print "\n*** When on source, adjust IF powers"
                if useOsf and band in ['327', '430']:
                    print "Set IF2 Auto-adjust offset to 18 for real-time cyclic spectroscopy hardware."
                
                if raw_input("\n*** Press enter when on source and IF powers are adjusted. (Type 'skip' to skip adjusting puppi attenuators)") != 'skip':
                    logging.info("auto leveling puppi")
                    while True:
                        autoPuppiAtten()
                        if raw_input("\n*** Check guppi_adc_hist. If happy with levels, press enter. Type anything else to retry setting attenuators.") == '':
                            break
    
            par = psr_par(parfile)
            try:
                psr = par.PSRJ
            except:
                try:
                    psr = par.PSR
                except:
                    raise Exception( "Error reading source name from parfile.")
            try:
                ra = par.RAJ
                dec = par.DECJ
            except:
                try:
                    ra = par.RA
                    dec = par.DEC
                except:
                    raise Exception( "Error reading RA/DEC from parfile.")
                
            obsmode = "COHERENT_CAL"
            tscan = 90.0
            polarg = ""
            acc_len = 1
            rcvrconfig = rcvrConfigs[band]
            nchan = rcvrconfig['nchan']
            bw = rcvrconfig['bw']
            freq = rcvrconfig['freq']
            pol = rcvrconfig['pol']
            rcvr = rcvrconfig['rcvr']
            
            if lastBand != band:
                if band == "327":
                    print "\n\n*** This is the 327 MHz receiver, so issue 'setcal hcorcal' and 'cal25 on'\n"
                else:
                    print "\n\n*** Issue 'setcal lcorcal' and 'cal25 on'\n"
            if raw_input("Press enter when finished or type 'skip' to skip the cal scan") != 'skip':
                
                logging.info("Configuring puppi params for cal")
                        
                # Build guppi_set_params call and execute it
                cmd = "guppi_set_params -f -i -b 2048 -t 10.0 -m %s %s" % (obsmode,polarg)
                cmd += " --nchan=%d --bw=%.1f --freq=%.1f -T %.1f" % (nchan, bw, freq, tscan)
                cmd += " --acc_len=%d --projid=%s --datadir=/data/puppi/%s" % (acc_len,projid,projid)
                cmd += " --src=%s --ra=%s --dec=%s --feed_pol=%s --frontend=%s" % (psr,ra,dec,pol,rcvr)
                cmd += " -P %s" % parfile
                logging.debug(cmd)
                
                os.system(cmd)
                
                if useOsf:
                    osf.stopData()
                    osf.killAll(osfnodes)
                logging.info("Starting calscan on %s" % (psr))
                time.sleep(1)
                gup.send("GPUS/DAQ/server","START")
                time.sleep(7)
                gup.arm()
                calend = time.time() + 90.0
                logging.debug("Puppi armed")
                
                print "Waiting for cal scan to finish",
                while gup.get("GPU5/DAQ/DAQSTATE") != 'stopped':
                    print ".",
                    sys.stdout.flush()
                    time.sleep(5)
                
                print "\n"
                logging.info("cal scan finished")
            lastBand = band

            
    if test:
        timeStart = startAt
    else:
        timeStart = time.time()
    
    riseSetList,startTimes,stopTimes,scanEndTimes = observe_utils.generateObservingPlan(session,timeStart,endAt,minTimePerSource=minTimePerSource)    
    
            
    print "\n=== Observing Plan ===\n"
    for n,(source,parfile,band) in enumerate(session):
        if n == 0:
            tstart = timeStart
             
        print "Observing",source,"at",band,"from",time.ctime(tstart),"to",time.ctime(scanEndTimes[n]), ("(%.1f mins)" % ((scanEndTimes[n]-tstart)/60.0))
        tstart = scanEndTimes[n]
        
    print "\n"
    print timeStart,time.time()
    while timeStart > time.time() and not test:
        print ("\r%.1f minutes until start time" % ((timeStart -time.time())/60.0)),
        sys.stdout.flush()
        time.sleep(10)
    print "\nstarting..."
    for (source,parfile,band),endTime in zip(session,scanEndTimes):
        logging.info("\n\n*** Slew to source: "+ source)
        sourceStartTime = time.time()
        timeForThisSource = endTime - sourceStartTime
        print "\n*** Select band:", band
        if raw_input("Press enter when finished. (Type 'skip' to skip configuring puppi)") != 'skip':
            mode,clk = guppiConfigs[band]
            if mode != lastMode:
                logging.info("Configuring puppi with mode " + mode +" and clock " + clk)
                gup.mode(mode)
                gup.reset(clk)
                gup.arm()
                lastMode = mode
            else:
                logging.info("previous PUPPI mode was the same as current, so not reconfiguring. Mode is: " + mode)

        print "\n*** When on source, adjust IF powers (auto-level IF1, then auto-level IF2)"
        if useOsf and band in ['327', '430']:
            print "Set IF2 Auto-adjust offset to 18 for real-time cyclic spectroscopy hardware."
        
        if raw_input("\n*** Press enter when on source and IF powers are adjusted. (Type 'skip' to skip adjusting puppi attenuators)") != 'skip':
            logging.info("auto leveling puppi")
            while True:
                autoPuppiAtten()
                if raw_input("\n*** Check guppi_adc_hist. If happy with levels, press enter. Type anything else to retry setting attenuators.") == '':
                    break

        par = psr_par(parfile)
        try:
            psr = par.PSRJ
        except:
            try:
                psr = par.PSR
            except:
                raise Exception( "Error reading source name from parfile.")
        try:
            ra = par.RAJ
            dec = par.DECJ
        except:
            try:
                ra = par.RA
                dec = par.DEC
            except:
                raise Exception( "Error reading RA/DEC from parfile.")
            
        obsmode = "COHERENT_CAL"
        tscan = 90.0
        polarg = ""
        acc_len = 1
        rcvrconfig = rcvrConfigs[band]
        nchan = rcvrconfig['nchan']
        bw = rcvrconfig['bw']
        freq = rcvrconfig['freq']
        pol = rcvrconfig['pol']
        rcvr = rcvrconfig['rcvr']
        
        if band == "327":
            print "\n\n*** This is the 327 MHz receiver, so issue 'setcal hcorcal' and 'cal25 on'\n"
        else:
            print "\n\n*** Issue 'setcal lcorcal' and 'cal25 on'\n"
        if raw_input("Press enter when finished or type 'skip' to skip the cal scan") != 'skip':
            
            logging.info("Configuring puppi params for cal")
            
    
            
            # Build guppi_set_params call and execute it
            cmd = "guppi_set_params -f -i -b 2048 -t 10.0 -m %s %s" % (obsmode,polarg)
            cmd += " --nchan=%d --bw=%.1f --freq=%.1f -T %.1f" % (nchan, bw, freq, tscan)
            cmd += " --acc_len=%d --projid=%s --datadir=/data/puppi/%s" % (acc_len,projid,projid)
            cmd += " --src=%s --ra=%s --dec=%s --feed_pol=%s --frontend=%s" % (psr,ra,dec,pol,rcvr)
            cmd += " -P %s" % parfile
            logging.debug(cmd)
            
            
            os.system(cmd)
        
            if band not in ['327','430'] or not useOsf:
                if useOsf:
                    osf.stopData()
                    osf.killAll(osfnodes)
                logging.info("Starting calscan on %s" % (psr))
                time.sleep(1)
                gup.send("GPUS/DAQ/server","START")
                time.sleep(7)
                gup.arm()
                calend = time.time() + 90.0
                logging.debug("Puppi armed")
            else:
                logging.info("Starting calscan on %s with OSF+Puppi" % (psr))
                # start puppi
                time.sleep(1)
                servers = [('GPU%d/DAQ/server' %x) for x in puppinodes]
                cmds = ['START' for x in puppinodes]
                gup.send(servers,cmds)
                time.sleep(7)
                gup.arm()
                try:
                    scannum = int(gup.get("GPU5/DAQ/SCANNUM"))
                except:
                    logging.warn("Couldn't get scannumber from guppi status, so just auto incrementing cyclic spectroscopy scan number")
                    scannum = None
                calend = time.time() + 90.0
                logging.debug("Puppi armed")
                osf.setupObs(parfile, int(band), gpus=osfnodes, projid=projid)
                osf.killAll(osfnodes)
                osf.monitor(osfnodes)
                osf.autoAtten()
                try:
                    osf.autoLevel(gpus=osfnodes, pass2=False)
                except Exception, e:
                    logging.warn("Per-channel auto leveling failed, proceeding anyway. Error was %s" % str(e))
                tdspsr = calend - time.time()
                if tdspsr < 10:
                    logging.warning("Ran out of time for dspsr cal")
                else:
                    logging.info("starting %.1f minute dspsr cal scan" % (tdspsr/60.0))
                    osf.startDspsr(gpus=osfnodes, cal = True, scanlen=tdspsr,ncyclic=256,nbin=256,extra='-nsub 5',lsubint=10.0, countdown=False, scannum = scannum)
    
            print "Waiting for cal scan to finish",
            while gup.get("GPU5/DAQ/DAQSTATE") != 'stopped':
                print ".",
                sys.stdout.flush()
                time.sleep(5)
            
            print "\n"
            logging.info("cal scan finished")
        
        print "\n\n*** Turn OFF cal: issue 'cal25 off'"
        if raw_input("Press enter when finished (type 'skip' to skip the pulsar scan)") != 'skip':
            
            ################# On to real scan
            
            tscan = endTime - time.time() - 10.0 # deduct 10 seconds for arming etc.
            if tscan < 0:
                raise Exception("Error! ran out of time! scan time is negative!")
            obsmode = "COHERENT_FOLD"
            
            # Build guppi_set_params call and execute it
            cmd = "guppi_set_params -f -i -b 2048 -t 10.0 -m %s %s" % (obsmode,polarg)
            cmd += " --nchan=%d --bw=%.1f --freq=%.1f -T %.1f" % (nchan, bw, freq, tscan)
            cmd += " --acc_len=%d --projid=%s --datadir=/data/puppi/%s" % (acc_len,projid,projid)
            cmd += " --src=%s --ra=%s --dec=%s --feed_pol=%s --frontend=%s" % (psr,ra,dec,pol,rcvr)
            cmd += " -P %s" % parfile
            logging.debug(cmd)
            
            os.system(cmd)
            
            if band not in ['327','430'] or not useOsf:
                if useOsf:
                    osf.stopData()
                    osf.killAll(osfnodes)
                logging.info("Starting %.1f minute observation on %s" % (tscan/60.0, psr))
                time.sleep(1)
                gup.send("GPUS/DAQ/server","START")
                time.sleep(7)
                gup.arm()
                logging.debug("Puppi armed")
            else:
                logging.info("Starting %.1f minute observation on %s with OSF+Puppi" % (tscan/60.0, psr))
                # start puppi
                time.sleep(1)
                servers = [('GPU%d/DAQ/server' %x) for x in puppinodes]
                cmds = ['START' for x in puppinodes]
                gup.send(servers,cmds)
                time.sleep(7)
                gup.arm()
                logging.debug("Puppi armed")
                try:
                    scannum = int(gup.get("GPU5/DAQ/SCANNUM"))
                except:
                    logging.warn("Couldn't get scannumber from guppi status, so just auto incrementing cyclic spectroscopy scan number")
                    scannum = None
                osf.setupObs(parfile, int(band), gpus=osfnodes, projid=projid)

                tdspsr = endTime - time.time()
                logging.info("starting %.1f minute dspsr scan" % (tdspsr/60.0))
                osf.startDspsr(gpus=osfnodes, scanlen=tdspsr,ncyclic=256,nbin=256,extra='-nsub 12',lsubint=10.0,countdown=False,scannum=scannum)
            
            logging.info("%.1f minute observation on %s should finish around %s" % (tscan/60.0, psr, time.ctime(time.time()+tscan)))
            print "waiting for scan to finish"
            
            while gup.get("GPU5/DAQ/DAQSTATE") != 'stopped':
                print ".",
                sys.stdout.flush()
                time.sleep(5)
            
            print "\n"
            logging.info("data scan finished")
            
    if useOsf:
        logging.info("Finishing up; making sure cyclic spectrometer is turned off.")
        osf.stopData()
        osf.killAll(osfnodes)
        

### Puppi Atten stuff (break out into new module)

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
    
    
    
def epochToMJD(t=None):
    """
        Return the current MJD accurate to ~1 sec if no argument
        otherwise convert argument from unix epoch to MJD
    """
    from pyslalib import slalib as s
    if t is None:
        YY, MM, DD, hh, mm, ss, wday, yday, isdst = time.gmtime()
    else:
        YY, MM, DD, hh, mm, ss, wday, yday, isdst = time.gmtime(t)
    mjd_f, J = s.sla_dtf2d(hh, mm, ss)
    mjd_i, J = s.sla_cldj(YY, MM, DD)
    return  mjd_i + mjd_f


    
if __name__ == "__main__":
    calSession = []
    test = False
    minTimePerSource = 2200.0
    print sys.argv
    if len(sys.argv) != 2:
        print "please specify session definition file"
        sys.exit(1)
    execfile(sys.argv[1])
    observe(session,endAt,projid,calSession=calSession,test=test)
    
    