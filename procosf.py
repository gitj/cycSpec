import numpy as np
import glob
import os
import scipy.signal
import scipy.optimize
from matplotlib import pyplot as plt
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
from matplotlib.figure import Figure                         
from matplotlib.backends.backend_agg import FigureCanvasAgg 


def interlace(x1,x2,ncyc=512):
    xx = np.empty_like(x1)
    nblock = xx.shape[0]/ncyc
    chk = ncyc/2
    offs = ncyc/4
    xx[:,:] = x1[:,:]
    for b in range(1,nblock-1):
        start = offs + ncyc*b
        stop = start + chk
        print b, start, stop, " : ", start+chk, stop+chk
        xx[(start+chk):(stop+chk),:] = x2[start:stop,:]
    return xx

def renameBrokenFiles():
    raws = glob.glob('/data/puppi/P2721/gpu*/cspuppi_*.raw')
    fixed = set(glob.glob('/data/puppi/P2721/gpu*/cspuppi_*.rawf'))
    torename = [x for x in raws if (x+'f') in fixed]
    for fn in torename:
        os.rename(fn,fn+'.broken')
    return torename
    

def getAllSubints(fglob='/data/puppi/P2721/gpu*/cspuppi_*_0097.0000.*.ar.fix',subints=(0,3),band=327):
    fnames = glob.glob(fglob)
    if not fnames:
        raise Exception("no files")
    print "loading:",fnames
    return getSubints(fnames, subints=subints,band=band)

def getSubints(fnames,subints=(0,3),band=327):
    import psrchive
    cfsegs = dict([(getObsFreqs(s, band=band, x16=True, Fs=1024.0)[16],s) for s in range(8)])
    ar0 = psrchive.Archive_load(fnames[0])
    npfb = 32
    nch = ar0.get_nchan()
    nbins = ar0.get_nbin()
    npol = ar0.get_npol()
    source = ar0.get_source()
    nsubints = subints[1]-subints[0] 
    i0 = ar0.get_Integration(subints[0])
    epoch = i0.get_epoch().in_days()
    ref_freq = 1/i0.get_folding_period()
    del ar0
    ncyc = nch/npfb
    data = np.zeros((nsubints,npol,4*nch + ncyc/2,nbins), dtype='float32')
    frqs = np.zeros((4*nch + ncyc/2,))
    for fname in fnames:
        ar = psrchive.Archive_load(fname)
        d = ar.get_data()
        s = cfsegs[ar.get_centre_frequency()]
        chk = s/2
        if chk < 2:
            chk += 2
        else:
            chk -= 2
        half = s % 2
        for ch in range(npfb):
            s1 = ncyc/4 + ch*ncyc
            s2 = s1 + ncyc/2
            d1 = s1 + half*ncyc/2 + ncyc/4 + chk*nch
            d2 = d1 + ncyc/2     
            print "s",s,"chk",chk,"s1:2",s1,s2,"d1:2",d1,d2
            if s1 < d.shape[2] and s2 < d.shape[2]:
                try:
                    data[:,:,d1:d2,:] = d[subints[0]:subints[1],:,s1:s2,:]
                except Exception, e:
                    print "problem with", fname, s
                    print e
            for ss,dd in zip(range(s1,s2),range(d1,d2)):
                if ss >= nch:
                    continue
                if dd >= frqs.shape[0]:
                    continue
                frqs[dd] = ar.get_Integration(0).get_Profile(0,ss).get_centre_frequency()
        del ar
    valid = np.flatnonzero(frqs)
    frqs = frqs[valid[0]:valid[-1]+1]
    return data[:,:,valid[0]:valid[-1]+1,:], frqs, epoch, ref_freq, source

def cleanAndNormalize(on,off,discard=0.01):
    #data dimensions are subint, poln, channel
    ds = on/off - 1.0
    if np.any(~np.isfinite(ds[-1,:])):
        ds = ds[:-1,:]
    shape = ds.shape
    ndat = np.prod(shape)
    ds.shape = (ndat,)
    ordering = ds.argsort()
    ds[ordering[:int(ndat*discard)]] = 0.0
    ds[ordering[-int(ndat*discard):]] = 0.0
    ds.shape = shape
    return ds.mean(1)  #scrunch the polns

pulseRegions = {('J2229+2643',327): ((100,180),(200,250)),
                ('J2229+2643',430): ((0,45),(100,200)),
                ('B1953+29',430): ((176,207),(50,100)),
                ('B1953+29',1400): ((207,222),(75,150)),
                ('2317+1439',327): ((105,130),(200,250)),
                ('2317+1439',430): ((111,140),(200,250)),
                ('1944+0907',430): ((187,255),(100,150)),
                ('J2043+1711',430): ((128,156),(200,250)),
                ('2017+0603', 430): ((220,255),(20,100)),
                ('1853+1303', 430): ((150,200),(40,120)),
                ('0023+0923', 430): ((51,75),(100,250)),
                ('0030+0451', 430): ((0,36), (175,225)),
                ('1640+2224', 430): ((223,250), (50,200)),
                ('1855+09', 430): ((212,246), (10,90)),
                }

def doAllPlots():
    pngs = {}
    for sn in range(113,156):
        try:
            ds,acf,frqs,times,source,band,epoch,fname = plotDynSpecForScan(sn)
            if not pngs.has_key(source):
                pngs[source] = {}
            if not pngs[source].has_key(band):
                pngs[source][band] = []
            pngs[source][band].append((epoch,fname))
        except Exception,e:
            print "Error on scan",sn
            print str(e)
    for source in pngs.keys():
        sourcepngs = []
        bands = pngs[source].keys()
        bands.sort()
        for band in bands:
            fns = pngs[source][band]
            epochsorted = sorted(fns, key = lambda x: x[0])
            sourcepngs = sourcepngs + [x[1] for x in epochsorted]
            
        makePdf(sourcepngs,('/home/gjones/plots/%s.pdf' % source))
        

def makePdf(pngs,pdfname):
    pngstr = ' '.join(pngs)
    cmd = "convert %s book.mng" % pngstr
    print cmd
    os.system(cmd)
    cmd = "convert book.mng -compress jpeg -quality 50 %s" % pdfname
    print cmd
    os.system(cmd)

def plotDynSpecForScan(sn,mjd):
    import psrchive
    files = glob.glob('/lakitu/data/P2721/gpu*/rt/rtcs_%d*_%04d_*%04d.ar.fix' % (mjd,sn,1))
    if len(files) == 0:
        print "found no files for scan",sn
        return
    ar = psrchive.Archive_load(files[0])
    source = ar.get_source()
    epoch = ar.get_Integration(0).get_start_time().in_days()
    cf = ar.get_centre_frequency()
    del ar
    if cf < 370:
        band = 327
    elif cf > 370 and cf < 500:
        band = 430
    else:
        band = 1400    
    try:
        onp,offp = pulseRegions[(source,band)]
    except:
        print "Could not find pulse regions for source", source,"at band", band
        return
    print "found source", source,"band", band, "for file", files[0]
    print "getting dynamic spectrum data"
    on,off,frqs,times = getDynSpec(('/lakitu/data/P2721/gpu*/rt/rtcs_%d*_%04d_' % (mjd,sn))+ '*%04d.ar.fix', band, onp, offp)
    
    # XXX Hack to avoid bug in ACF code for odd number of subints
    if on.shape[0] % 2 == 1:
        on = on[:-1]
        off = off[:-1]
        times = times[:-1]
    print "cleaning and normalizing"
    ds = cleanAndNormalize(on, off)
    times= times[:ds.shape[0]]  #remove last time entry if cleaner killed last subint
    print "computing ACF"
    acf = computeAcf(ds)
    fig = Figure(figsize=(10,12))
    fig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
    plotDynSpecAcf(ds,acf,frqs,times,fig=fig)
    fname = os.path.join('/home/gjones/plots',(os.path.split(files[0])[1] + '_'+ 'dynspec.png'))
    print fname
    esc_fname = fname.replace('_',r'\_')
    fig.suptitle(('%s @ %d MHz AO %s' % (source,band,esc_fname)),size='medium')
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(fname)
    return ds,acf,frqs,times,source,band,epoch,fname

def plotDynSpecAcf(ds,acf,frqs,times,fig = None):
    times = times - times[0]
    df = np.abs(frqs[1]-frqs[0])
    dt = np.diff(times).mean()
    pt,t,tacf,pf,fr,facf,offset1,scale1 = fitAcf(acf, dt=times[1], df=df)
    try:
        pt2,t2,tacf2,pf2,fr2,facf2,offset2,scale2 = fitAcf(acf, dt=times[1], df=df,maxf = 1.0)
    except:
        pf2 = (1.0,)
        offset2 = 0
        scale2 = 0.1
    fscale = int(5*pf[0]/df)
    print fscale
    if fscale < 512:
        fscale = 512
    if fig is None:
        f = plt.figure()
    else:
        f = fig
    ax1 = f.add_subplot(3,2,1)
    ax2 = f.add_subplot(3,2,2)
    ax3 = f.add_subplot(3,1,2)
    ax4 = f.add_subplot(3,1,3)
    ax1.plot(t,tacf)
    ax1.plot(t,gammaIs(t,pt[0]), linewidth=2, alpha=0.5, label=((r'$\exp{(-(t/%.2f)^{5/3})}$' % pt[0])))
    ax1.legend()
    ax1.set_title('Normalized temporal ACF')
    ax1.set_xlabel('Lag (s)')
    ax2.plot(fr,facf)
    ax2.plot(fr,gammaIf(fr,pf[0]), linewidth=2, alpha=0.5, label=((r'$\exp{(\frac{-(\log{2})\Delta f}{%.3f})}$' % pf[0])))
    ax2.plot(fr,(gammaIf(fr,pf2[0])*scale2 + offset2-offset1)/scale1, linewidth=2, alpha=0.5, label=((r'$\exp{(\frac{-(\log{2})\Delta f}{%.3f})}$' % pf2[0])))
    ax2.set_xlim(0,fscale*df)
    ax2.legend()
    ax2.set_title('Normalized frequency ACF')
    ax2.set_xlabel('Lag (MHz)')
    ds2 = rebinAxis(ds, 2048, axis=1)
    ax4.imshow(ds2,aspect='auto',origin='lower',extent=[frqs[0],frqs[-1],times[0],times[-1]])
    ax4.text(0.1,0.9,"Dynamic Spectrum",transform=ax4.transAxes,bbox=dict(facecolor='w',alpha=0.5))
    ax4.set_xlabel('Freq (MHz)')
    ax4.set_ylabel('Time (s)')
    
    
    if fscale > 1024:
        scaleby = int(np.ceil(fscale/1024))
        acf2 = rebinAxisBy(acf, scaleby, axis=1)
    else:
        acf2 = acf
        scaleby = 1
    print fscale,scaleby
    acf2 = acf2/acf2[tuple(np.array(acf2.shape)/2 + 5)]
    acf2 = acf2 - acf2[int(0.75*acf2.shape[0]):,int(0.75*acf2.shape[1]):].mean()
#    acfs = acf2.flatten()
#    acfs.sort()
#    levels = [acfs[int(acfs.shape[0]*x)] for x in (1-10**-np.linspace(2,5,10))]
#    levels = np.linspace(0.4,1,6)
    nf = acf2.shape[1]    
    extent = [-df*fscale, df*fscale, -dt*acf2.shape[0]/2.0,dt*acf2.shape[0]/2.0]
    im = ax3.imshow(acf2[:,nf/2-fscale/scaleby:nf/2+fscale/scaleby],aspect='auto',origin='lower',extent=extent)
    
    smoothfactor =8
    acfsm = rebinAxisBy(acf2,smoothfactor,axis=1)
    levels = np.array([.5,.75,.9,.95,.99])
    levels = np.linspace(.1,1,10)
    nfsm = acfsm.shape[1]
    ax3.text(0.1,0.9,"2-D ACF",transform=ax3.transAxes,bbox=dict(facecolor='w',alpha=0.5))
    im.set_clim(-.1,1)
    fextent = fscale/(scaleby*smoothfactor)
    cs = ax3.contour(acfsm[:,nfsm/2-fextent:nfsm/2+fextent],extent=extent,
                     levels = levels, cmap=plt.cm.gray_r, linewidths=1)
    cb = f.colorbar(im,ax=ax3)
    try:
        cb.add_lines(cs)
    except ValueError:
        pass
    
    ax3.set_xlabel('Lag (MHz)')
    ax3.set_ylabel('Lag (s)')

def computeAcf(in1,correct=True):
    s1 = np.array(in1.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex))
    size = 2**np.ceil(np.log2(s1*2))
    IN1 = np.fft.fftn(in1,size)
    IN1 *= np.conj(IN1)
    ret = np.fft.ifftn(IN1)
    del IN1
    if not complex_result:
        ret = ret.real
    osize = s1
    output = scipy.signal.signaltools._centered(np.fft.fftshift(ret),osize)
    if not correct:
        return output
    corrections = []
    for d in range(output.ndim):
        ndat = output.shape[d]
        c = np.zeros((ndat,))
        half = ndat % 2
        c[:ndat/2] = (ndat*1.0/(ndat-np.arange(1,ndat/2+1)))[::-1]
        c[ndat/2:] = ndat*1.0/(ndat-np.arange(ndat/2 + half))
        corrections.append(c)
    return output * np.outer(corrections[0],corrections[1])
    
def gammaIs(x,a):
    return np.exp(-(x/a)**(5/3.0))

def objGammaIs(p,x,y):
    return (np.abs(y-gammaIs(x,*p))**2).sum()

def gammaIf(x,a):        
    return np.exp(-np.abs((np.log(2)*x/a)))

def objGammaIf(p,x,y):
    return (np.abs(y-gammaIf(x,*p))**2).sum()

def fitAcf(acf,dt=1.0,df=1.0,maxf=None):
    """
    assumes time axis first
    """
    tacf = acf[acf.shape[0]/2+1:,acf.shape[1]/2]
    tacf = tacf - tacf[tacf.shape[0]/2:].mean()
    tacf = tacf/tacf[0]
    times = dt * np.arange(1,tacf.shape[0]+1)
    x0 = times[np.abs(tacf-tacf.max()/2.0).argmin()]
    pt = scipy.optimize.fmin(objGammaIs,x0=[x0],args=(times,tacf))
    
    facf = acf[acf.shape[0]/2,acf.shape[1]/2+2:]
    freqs = df * np.arange(2,facf.shape[0]+2)
    if maxf:
        nf = int(maxf/df)
        facf = facf[:nf]
        freqs = freqs[:nf]
    offset = facf[facf.shape[0]/2:].mean()
    facf = facf - offset
    scale = facf[0]
    facf = facf/scale
    x0 = freqs[np.abs(facf-facf.max()/2.0).argmin()]
    pf = scipy.optimize.fmin(objGammaIf,x0=[x0],args=(freqs,facf))
    pf[0] = np.abs(pf[0])
    return pt,times,tacf,pf,freqs,facf,offset,scale
    

def computeCcf(in1,in2):
    """Correlate two N-dimensional arrays using FFT. See convolve.

    """
    s1 = np.array(in1.shape)
    s2 = np.array(in2.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex) or
                      np.issubdtype(in2.dtype, np.complex))
    size = s1+s2
    IN1 = np.fft.fftn(in1,size)
    IN1 *= np.conj(np.fft.fftn(in2,size))
    ret = np.fft.ifftn(IN1)
    del IN1
    if not complex_result:
        ret = ret.real
    if np.product(s1,axis=0) > np.product(s2,axis=0):
        osize = s1
    else:
        osize = s2
    return scipy.signal.signaltools._centered(ret,osize)

def rebinAxisBy(data,nfact,axis=0,oper=np.mean):
    nout = int(data.shape[axis]/nfact)
    return rebinAxis(data,nout,axis=axis,oper=oper)
    
def rebinAxis(data,nout,axis=0, oper = np.mean):
    nout = int(nout)
    nstart = data.shape[axis]
    nchk = nstart/nout
    nstop = nchk*nout
    index = [slice(None) for x in range(axis)] + [slice(nstop)] 
    if data.ndim > axis:
        index = index + [slice(None) for x in range(data.ndim-axis-1)]
    index = tuple(index)
    newshape = list(data.shape)
    newshape[axis] = nout
    newshape.insert(axis+1,-1)
    return np.mean(data[index].reshape(tuple(newshape)),axis=axis+1)

validFrequencies = {327: (295,370),
                    430: (420,450),
                    1400: (1000,1475)}
def getDynSpec(fglob='/data/puppi/P2721/gpu*/cspuppi_*_0097.%04d.*.ar.fix',band=327, onp=(100,150),offp=(200,250), discard=True):
    import psrchive
    fn = 0
    if len(glob.glob(fglob % fn)) == 0:
        fn = 1 
    blockson = []
    blocksoff = []
    frqblks = []
    times = []
    cfsegs = dict([(getObsFreqs(s, band=band, x16=True, Fs=1024.0)[16],s) for s in range(8)])
    nsubtot = 0
    while True:
        fnames = glob.glob(fglob % fn)
        if not fnames:
            break
        print "loading:",fnames
        ar0 = psrchive.Archive_load(fnames[0])
        npfb = 32
        nch = ar0.get_nchan()
        nsubint = ar0.get_nsubint()
        for k in range(nsubint):
            times.append(ar0.get_Integration(k).get_epoch().in_seconds())
            
        nsubtot += nsubint
        del ar0
        ncyc = nch/npfb
        dataon = np.zeros((nsubint,2,4*nch + ncyc/2), dtype='float32')
        dataoff = np.zeros((nsubint,2,4*nch + ncyc/2), dtype='float32')
        frqs = np.zeros((4*nch + ncyc/2,))
        for fname in fnames:
            ar = psrchive.Archive_load(fname)
            d = ar.get_data()
            don = d[:nsubint,:,:,onp[0]:onp[1]].mean(3)
            doff = d[:nsubint,:,:,offp[0]:offp[1]].mean(3)
            s = cfsegs[ar.get_centre_frequency()]
            chk = s/2
            if chk < 2:
                chk += 2
            else:
                chk -= 2
            half = s % 2
            for ch in range(npfb):
                s1 = ncyc/4 + ch*ncyc
                s2 = s1 + ncyc/2
                d1 = s1 + half*ncyc/2 + ncyc/4 + chk*nch
                d2 = d1 + ncyc/2     
#                print "s",s,"chk",chk,"s1:2",s1,s2,"d1:2",d1,d2
                if s1 < d.shape[2] and s2 < d.shape[2]:
                    try:
                        this_nsubint = d.shape[0]
                        dataon[:this_nsubint,:,d1:d2] = don[:,:,s1:s2]
                        dataoff[:this_nsubint,:,d1:d2] = doff[:,:,s1:s2]
                    except Exception, e:
                        print "problem with", fname, s
                        print e
                for ss,dd in zip(range(s1,s2),range(d1,d2)):
                    if ss >= nch:
                        continue
                    if dd >= frqs.shape[0]:
                        continue
                    frqs[dd] = ar.get_Integration(0).get_Profile(0,ss).get_centre_frequency()
            del ar
        blockson.append(dataon)
        blocksoff.append(dataoff)
        frqblks.append(frqs)
        fn += 1
    on = np.concatenate(blockson,axis=0)[:,:,256:]
    off = np.concatenate(blocksoff,axis=0)[:,:,256:]
    frqs = frqblks[0][256:]
    nonz = np.flatnonzero(off.min(0).min(0))
    start = nonz[0]
    stop = nonz[-1]+1
    on = on[:,:,start:stop]
    off = off[:,:,start:stop]
    frqs = frqs[start:stop]
    if discard:
        validf = validFrequencies[band]
        print validf
        start = np.flatnonzero(frqs>=validf[0])[0]
        stop = np.argmin(np.abs(frqs- validf[1]))
#        stop = np.flatnonzero(frqs<= validf[1])[-1]
        print start,stop
        on = on[:,:,start:stop]
        off = off[:,:,start:stop]
        frqs = frqs[start:stop]
        
        
    # find longest contiguous frequency chunk
    difs = np.diff(np.hstack(([0],np.diff(frqs),[0])))
    starts, = np.where(difs>0)
    ends, = np.where(difs<0)
    runs = ends-starts
    start = starts[runs.argmax()]
    stop = ends[runs.argmax()]+1
    on = on[:,:,start:stop]
    off = off[:,:,start:stop]
    frqs = frqs[start:stop]
    
    times = np.array(times)
    return on,off,frqs,times
        

def fixApril12Files(fglob='/data/gpu/partial/gpu*/cspuppi/*.raw'):
    files = glob.glob(fglob)
    errors = {}
    for fn in files:
        print "fixing", fn
        try:
            fixRawFile(fn)
        except Exception, e:
            print "failed to fix", fn
            print "error was:",e
            errors[fn] = e
    return errors

def dspsrCal(fglob='/data/gpu/partial/gpu*/cspuppi/cs*.raw'):
    fnames = glob.glob(fglob)
    for fname in fnames:
        print fname
        fout = fname + '.cal'
        hdr = getRawHeader(fname)
        cal = hdr['CAL_MODE']
        if cal.find('ON') >= 0:
            cmd = 'source /home/gpu/gjones/puppi.sh; dspsr -overlap -O %s -c0.04 -D0.00000001 -j "e type=PolnCal" -x65536 -cuda 0,0 -F32:D -L10.0 -b256 -A %s &> %s.log' % (fout,fname,fout)
            print cmd
            os.system(cmd)
            
def rmScans(scannos=range(75,92)):
    for sn in scannos:
        fnames = (glob.glob('/data/gpu/partial/gpu*/cspuppi/rt/rt*_%04d_*' % sn) 
                  + glob.glob('/data/gpu/partial/gpu*/cspuppi/cspuppi*_%04d_*' % sn))
        for fname in fnames:
            print fname
            os.unlink(fname)
    
def dspsrGlob(fglob='/data/gpu/partial/gpu*/cspuppi/cs*.raw'):
    fnames = glob.glob(fglob)
    froots = list(set([x[:-10] for x in fnames]))
    failures = []
#    for froot in froots:
    for fname in fnames:
        print fname
        fout = fname + '.cyc'
        cmd = 'source /home/gpu/gjones/puppi.sh; dspsr -cont -overlap -O %s -E /home/gpu/tzpar/B1937+21.par -x65536 -cuda 0,0 -cyclic 512 -F32:D -L10.0 -b256 -d2 -A %s &> %s.log' % (fout,fname,fout)
        print cmd
        os.system(cmd)
        try:
            fixCyclicDedisp(fout+'.ar')
        except:
            failures.append(fname)
    return failures

def fixRawFile(fn):
    hdr = getRawHeader(fn)
    seg = int(hdr['SEGMENT'])
    band = int(hdr['FRONTEND'][1:-1].strip()) #won't work for lwide
    cf = getObsFreqs(seg, band, x16=True, Fs=1024.0)[16]
    fixRawFreqs(fn,cf)
    
    
def fixCyclicDedisp(fname, nchan=32, overwrite=False, ext='fix'):
    # copied from paul's fix_cyclic_dedisp script
    import psrchive
    import os
    import psr_utils
    arch = psrchive.Archive_load(fname)
    cf = arch.get_centre_frequency()
    bw = arch.get_bandwidth()
    f_lo = cf - bw/2.0
    nch = arch.get_nchan()
    pfb_nch = nchan
    pfb_bw = bw / pfb_nch
    chan_per_pfb = nch / pfb_nch
    dm = arch.get_dispersion_measure()
    for isub in range(arch.get_nsubint()):
        sub = arch[isub]
        per = sub.get_folding_period()
        for ichan in range(nch):
            pfb_cf = f_lo + ((ichan/chan_per_pfb)+0.5)*pfb_bw
            dt = psr_utils.delay_from_DM(dm,pfb_cf) - psr_utils.delay_from_DM(dm,cf)
            for ipol in range(sub.get_npol()):
                prof = sub.get_Profile(ipol,ichan)
                prof.rotate_phase(dt/per)
    #arch.set_dedispersed(True) # doesn't work, lame...
    if (overwrite):
        outf = fname
    else:
        outf = fname + '.' + ext
    arch.unload(outf)
    os.system("psredit -m -c dmc=1 %s" % outf)
    

def getRawHeader(fn):
    fh = open(fn,'rb')
    lines = 0
    hdr = {}
    while lines < 100:
        x = fh.read(80)
        if x.startswith('END'):
            break
        key,value = x.split('=')
        hdr[key.strip()] = value.strip()
        lines += 1
    else:
        fh.close()
        raise Exception("Could not find segment in first 100 lines")
    fh.close()
    return hdr
    

def fixRawFreqs(fn,cf,bw=16.0,chbw=0.5):
    bs = 134223968
    fh = open(fn,'rb')
    fout = open(fn+'f','wb')
    fh.seek(0,2)
    nblk = fh.tell()/bs
    fh.seek(0)
    print "found %d blocks" % nblk
    for blk in range(nblk):
        print "block", blk
        fh.seek(blk*bs)
        data = fh.read(bs)
        offset = 0
        while True:
            x = data[offset*80:(offset+1)*80]
            if x.startswith('END'):
                print "found end at offset",offset
                fout.write(x)
                offset += 1
                break
            if x.startswith('OBSFREQ'):
                print "found"
                print x
                x = '%-8s=%21s' % ('OBSFREQ',str(cf))
                x = '%-80s' % x
                print "replacing with"
                print x
            elif x.startswith('CHAN_BW'):
                print "found"
                print x
                x = '%-8s=%21s' % ('CHAN_BW',str(chbw))
                x = '%-80s' % x
                print "replacing with"
                print x
            elif x.startswith('OBSBW'):
                print "found"
                print x
                x = '%-8s=%21s' % ('OBSBW',str(bw))
                x = '%-80s' % x
                print "replacing with"
                print x
            fout.write(x)
            offset += 1
        print "writing data block... ",
        fout.write(data[offset*80:])
        print "finished"
    fh.close()
    fout.close()

nyqdict = {327 : 2,
           430 : 2,
           1400 : 3,}
feparams = {327 : dict(FRONTEND='327', FD_POLN='LIN'),
            430 : dict(FRONTEND='430', FD_POLN='CIRC'),
            1400 : dict(FRONTEND='L-wide',FD_POLN='LIN')}
freqoffs = {327 : 1077,
            430 : 1180,
            1400 : 180}


def segFreq(s,Fs=1024.0,nyq=2,x16 = False):
    chwid = 2*Fs/8.0/512.0
    nch = 512
    if nyq == 2:
        fcs = Fs - (Fs/4.0 + np.fft.fftshift((chwid/2.0)*(np.arange(nch)-nch/2)) - chwid/2.0)
    elif nyq == 3:
        fcs = Fs + (Fs/4.0 + np.fft.fftshift((chwid/2.0)*(np.arange(nch)-nch/2)) - chwid/2.0) # check this one
    if not x16:
        #original x8 distributor case
        chunk = s/2
        half = s%2
        chans = range(chunk*nch/4, (chunk+1)*nch/4)[half::2]
        return fcs[chans], np.mod(np.array(chans,dtype='int')-1,512) #[chunk*nch/4:(chunk+1)*nch/4][half::2]
    else:
        if s <  4:
            #first 4 segments are not remapped
            half = (s%2)
            offset = 1
            pass
        elif s >=8:
            raise Exception("Non existant segments are not handled")
        else:
            # last 4 segments are remmaped to end
            s = s + 8   
            half = s%2
            offset = 1
        chunk = s/2
#        half = (s%2)
        chans = range(chunk*nch/8, (chunk+1)*nch/8)[half::2]
        return fcs[chans], np.mod(np.array(chans,dtype='int')-offset,512) #[chunk*nch/8:(chunk+1)*nch/8][half::2]

def getObsFreqs(seg,band,x16=True,Fs=1024.0):
    nyq = nyqdict[band]
    frqs,chans = segFreq(seg,Fs=Fs,nyq=nyq,x16=x16)
    if nyq == 2:
        frqs = freqoffs[band] - frqs
    else:
        frqs = freqoffs[band] + frqs
    return frqs
