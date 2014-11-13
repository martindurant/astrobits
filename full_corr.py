import general_rebin
import startup as s
from scipy.optimize import leastsq
import numpy as n

def gauss(x,x0,sigma):
    """gaussian which works on either scalar or array"""
    return exp( - ((x-x0) / (sqrt(2)*sigma))**2)

def corr(chunksize,chunkstep,time1,flux1,time2,flux2,timestep=0):
    """Do cross-correlation like general_rebin.chunkcrosscorr, but return
additional information: average fluxes, standard deviation and
auto-correlations"""
    mintimestep = timestep or max((n.diff(time1).min(),n.diff(time2).min()))
    mintime = max((time1.min(),time2.min()))
    maxtime = min((time1.max(),time2.max()))
    t = n.arange(mintime,maxtime,mintimestep)
    f1 = general_rebin.rebin(time1,flux1,t)
    f2 = general_rebin.rebin(time2,flux2,t)
    start = n.arange(mintime,maxtime-chunksize,chunkstep)
    nmax = int(chunksize/mintimestep)+1
    output = n.zeros((len(start),nmax))*1.
    output2 = n.zeros((len(start),nmax))*1.
    output3 = n.zeros((len(start),nmax))*1.
    ave1 = n.zeros(len(start))*1.
    ave2 = n.zeros(len(start))*1.
    std1 = n.zeros(len(start))*1.
    std2 = n.zeros(len(start))*1.
    for i in range(len(start)):
        ind = (t>start[i])*(t<=start[i]+chunksize)
        output[i,:ind.sum()] = general_rebin.correlate(f1[ind]-f1[ind].mean(),f2[ind]-f2[ind].mean(),'same')/((f1[ind].std())*(f2[ind].std())*len(f1[ind]))
        acf1 = general_rebin.correlate(f1[ind]-f1[ind].mean(),f1[ind]-f1[ind].mean(),'same')/(f1[ind].std()**2*len(f1[ind]))
        acf2 = general_rebin.correlate(f2[ind]-f2[ind].mean(),f2[ind]-f2[ind].mean(),'same')/(f2[ind].std()**2*len(f1[ind]))
        output2[i,:ind.sum()] = acf1
        output3[i,:ind.sum()] = acf2
        ave1[i] = f1[ind].mean()
        ave2[i] = f2[ind].mean()
        std1[i] = f1[ind].std()
        std2[i] = f2[ind].std()
    num = len(output[0])
    if num%2 == 0:
        lag = n.arange(-(num/2 - 1)*mintimestep,(num/2+0.01)*mintimestep,mintimestep)-mintimestep
    else:
        lag = n.arange(-(num-1)*mintimestep/2,(num-0.99)*mintimestep/2,mintimestep)
    return start,lag,output,output2,output3,ave1,ave2,std1,std2

def bart_err(acf1,acf2):
    """Derive empirical uncertainties on a CCF based on the two ACFs, following the Bartlett
formula (e.g., as nicely explaines in Smith & Vaughan, 2007). N is the number of pairs that
make each CCF point"""
    num = len(acf1)
    N = n.where(n.arange(1,num+1,1)<n.arange(num,0,-1),n.arange(1,num+1,1),n.arange(num,0,-1))
    sig2 = (acf1*acf2).sum()/N
    return n.sqrt(sig2)

def make_plot():
    logfile = "/home/durant/data/ultracam/ScoX1/18_run025_n_trm.red.log"
    ccd=1
    band="low"
    mjd,rel,err,mjd_utc,rate,e = ultraxte.correlate(logfile,ccd,root+logfiles[logfile]+mid+band+'.lc',fclip=2)
    start,lag,output,output2,output3,ave1,ave2,std1,std2=corr(50,30,mjd[rel>0.03]*24*3600,rel[rel>0.03],mjd_utc*24*3600,rate)
    imshow(output,extent=[lag[0]-0.26502320542931557/2,lag[-1]+0.26502320542931557/2,-15,start[-1]-start[0]+15])
    fitfunc = lambda p,x : abs(p[0]) * s.gauss(x,p[1],p[2])
    errfunc = lambda p,x,y : fitfunc(p,x) - y
    from fit import fit
    res=[]
    for i in range(len(start)):
        temp = fit(lag,output[i],bart_err(output2[i],output3[i]),p0[:],fitfunc,assume=0)
        res.append(temp)
    peak = n.ones(len(start))*1.0
    err = n.ones(len(start))*1.0
    try:
        peak[i] = res[i][0][1]
        err[i] = res[i][1][1]
    except:
        pass
    s.errorbar(peak2[peak2!=1],start[peak2!=1]-start[0],xerr=err2[peak2!=1],yerr=0,fmt='ko',capsize=0)
    s.xlabel("$\delta t$ (s)")
    s.ylabel("$T$ (s)")

fitfunc = lambda p,x : abs(p[0]) * gauss(x,p[1],p[2]) - abs(p[3]) * gauss(x,p[4],p[5])
errfunc = lambda p,x,y : fitfunc(p,x) - y

def fit_template(lag,template):
    ind1 = n.where((lag<0)*(template<0))[0][-1]
    ind0 = n.where(template[:ind1]>0)[0][-1]
    print lag[ind0],lag[ind1]
    ave1 = (-template*lag)[ind0:ind1+1].sum()*n.diff(lag).mean()
    amp1 = (-template)[ind0:ind1+1].sum()*n.diff(lag).mean()
    wid1 = (-template*(lag-ave1)**2)[ind0:ind1+1].sum()*n.diff(lag).mean()

    ind0 = n.where((lag>0)*(template>0))[0][0]
    ind1 = n.where(template[ind0:]<0)[0][0] + ind0
    print lag[ind0],lag[ind1]
    ave2 = (template*lag)[ind0:ind1].sum()*n.diff(lag).mean()
    amp2 = (template)[ind0:ind1+1].sum()*n.diff(lag).mean()
    wid2 = (template*(lag-ave2)**2)[ind0:ind1].sum()*n.diff(lag).mean()
    return ave1,amp1,wid1,ave2,amp2,wid2

def fit_ccf(lag,ccf,template):
    """Fit two-gaussian function to a single ccf (e.g., a single
line of corr). For use in anti-corr, plus-corr binary
optcial/X-ray CCFs"""
    cor = n.correlate(ccf,template,'same')
    ind = n.argmax(cor*(abs(lag)<15))
    return lag[ind],cor[ind]

def fit_acf(lag,acf):
    """Measure the width of an acf (i.e., a single line out
of corr."""
    ind = lag>0
    firstzero = argmax(acf[ind]<0)
    return acf[ind][0],(lag[ind][:firstzero]*acf[ind][:firstzero]).mean()
