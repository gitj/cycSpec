import numpy as np
import time
import os
try:
    from guppi_daq import astro_utils
except:
    import astro_utils
from parfile import psr_par

zaLimit = 19.6

def generateObservingPlan(session,startAt,endAt,minTimePerSource=2200.0):
    lastRise = startAt
    radecs = {}
    riseSet = {}
    riseSetList = []
    parfiles = []
    for source,parfile,band in session:
        if radecs.has_key(parfile):
            continue
        parfiles.append(parfile)
        par = psr_par(parfile)
        try:
            ra = par.RAJ
            dec = par.DECJ
        except:
            try:
                ra = par.RA
                dec = par.DEC
            except:
                raise Exception( "Error reading RA/DEC from parfile.")
            
        rh,rm,rs = map(float,ra.split(':'))
        dd,dm,ds = map(float,dec.split(':'))
        radeg = astro_utils.hms_to_rad(rh,rm,rs)*astro_utils.RADTODEG
        decdeg = dd + dm/60.0 + ds/3600.0
        radecs[parfile] = (radeg,decdeg)
        
        ts = np.arange(lastRise-2*3600,lastRise+5*3600,60.0) 
        zas = np.array([astro_utils.radec_to_azza(radeg,decdeg,epochToMJD(x),scope="ARECIBO")[1] for x in ts])
        #print zas[0],zas[-1]
        upTimes = ts[(zas<zaLimit)]
        keyHole = ts[(zas<2.0)]
        if len(keyHole):
            print "Warning:",parfile,"is in the key hole from",time.ctime(keyHole[0]),"to",time.ctime(keyHole[1])
        riseT,setT = (upTimes[0],upTimes[-1])
        riseSet[parfile] = riseT,setT
        lastRise = upTimes[0]
        riseSetList.append((riseT,setT))
        print "%10s: Rise %s  Set %s Up %.1f min." % (source[:10],time.ctime(riseT),time.ctime(setT), (setT-riseT)/60.)
    first = riseSetList[0]
    if first[0]< startAt:
        first = (startAt,first[1])
        begin = startAt
    else:
        begin = first[0]
    riseSetList[0] = first
    last = riseSetList[-1]
    if last[1]>endAt:
        last = (last[0],endAt)
    riseSetList[-1] = last
    ss,fv = optimizeSchedule(riseSetList, minTimePerSource = minTimePerSource)
    startTimes = ss[::2]
    stopTimes = ss[1::2]
    sourceEndTimes = dict(zip(parfiles,stopTimes))
    sourceStartTimes = dict(zip(parfiles,startTimes))
    nscans = {}
    for source,parfile,band in session:
        if nscans.has_key(parfile):
            nscans[parfile] += 1
        else:
            nscans[parfile] = 1
    
    scanEndTimes = []
    scanStart = begin
    lastpar = None
    tslew = 600.0
    for source,parfile,band in session:
        sourceEnd = sourceEndTimes[parfile]
        sourceStart = sourceStartTimes[parfile]
        if sourceStart < begin:
            sourceStart = begin
        timePerScan = (sourceEnd-sourceStart-tslew)/nscans[parfile]
        if lastpar == parfile:
            #we're on the same source
            scanEnd = scanStart + timePerScan
        else:
            #switching sources so add extra time for slewing
            scanEnd = sourceStart + timePerScan + tslew
        scanEndTimes.append(scanEnd)
        scanStart = scanEnd
        lastpar = parfile
    return riseSetList,startTimes,stopTimes,scanEndTimes

def optimizeSchedule(riseSetActual,minTimePerSource = 2200.,debug=False):
    from functools import partial
    t0 = riseSetActual[0][0]
    riseSet = [(x-t0,y-t0) for (x,y) in riseSetActual]
    nsources = len(riseSet)
    ntimes = nsources*2
    approxTimePerSource = (riseSetActual[-1][1] -riseSetActual[0][0])/nsources
    fvals = []
    def objective(x):
        rv =  -np.mean(x[1::2] - x[::2])
        avg = np.mean(x[1::2] - x[::2])
        rv = 0.05*np.abs((x[1::2] - x[::2])-avg).sum() + rv
        #rv = np.sum(np.abs((x[1::2] - x[::2])-approxTimePerSource))
        fvals.append(rv)
        return rv 
#        return 0 #-np.sum((np.diff(x)-approxTimePerSource)**2)
    
    x0 = np.zeros((ntimes,))
    x0[::2] = np.array([x[0]+.1 for x in riseSet])
    x0[1::2] = np.array([x[0]+minTimePerSource for x in riseSet])

    minTimeConstr = []
    riseConstr = []
    setConstr = []

    def cf(x,ns,mintime):
        #print ns
        rv =  (x[2*ns+1]-x[2*ns]) - mintime
        if debug and rv < -0.1:
            print "broke min time",(x[2*ns+1]-x[2*ns])
        return rv
    minTimeConstr = [partial(cf,ns=nsrc,mintime=minTimePerSource) for nsrc in range(nsources)]
    
    def causef(x,ns):
        rv = x[ns+1] - x[ns]
        if debug and rv < -0.1:
            print "broke causal constraint",ns, x[ns+1], x[ns]
        return rv
    causalConstr = [partial(causef,ns=nsrc) for nsrc in range(ntimes-1)]
    def risef(x,ns):
        #print ns
        rv = x[2*ns] - riseSet[ns][0]
        if debug and rv < -0.1:
            print "broke rise constr",ns,x[2*ns],riseSet[ns][0]
        return rv
            
    riseConstr= [partial(risef,ns=nsrc) for nsrc in range(nsources)]
    def setf(x,ns):
        #print ns
        rv = riseSet[ns][1] - x[2*ns+1]
        if debug and rv < -0.1:
            print "broke set constr",ns,x[2*ns+1],riseSet[ns][1]
        return rv
    setConstr= [partial(setf,ns=nsrc) for nsrc in range(nsources)]
    import scipy.optimize
    def cf2(x):
        rv = x[0]-riseSet[0][0]
        if debug and rv < -0.01:
            print "broke extra constr",rv
        return rv
    
    best = None
    mintime = minTimePerSource
    if debug:
        disp = 1
    else:
        disp = 0
    while True:
        minTimeConstr = [partial(cf,ns=nsrc,mintime=mintime) for nsrc in range(nsources)]

        x = scipy.optimize.fmin_cobyla(objective,x0,riseConstr + setConstr + minTimeConstr+causalConstr,maxfun=5000,rhobeg=10,disp=disp)#,rhoend=1e-2)#,rhobeg=2e-2,rhoend=1e-6)
        #print "checking constraints"
        cnstrs = ([cf2] + riseConstr + setConstr + minTimeConstr+causalConstr)
        bounds = np.array([f(x) for f in cnstrs])
        if np.any(bounds<-1):
            #some constraint failed
            break
        best = x
        mintime += 60.0
    if best is None:
        raise Exception("Failed to find a solution while optimizing the schedule!\nTry reducing the minTimePerSource parameter.")
    x = best
#    print x0
    x = x + t0
    print " "
    print "Total time observing:",np.sum(x[1::2]-x[::2])
    for ns in range(nsources):
        print "rise",time.ctime(riseSetActual[ns][0]),"start",time.ctime(x[2*ns])
        print "set",time.ctime(riseSetActual[ns][1]),"stop",time.ctime(x[2*ns+1])
        print "time for source",(x[2*ns+1]-x[2*ns])
    return x,np.array(fvals)
    

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
        
    