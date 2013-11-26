#from pylab import *
from Numeric import *
import Scientific.Statistics as stat
from Scientific.Functions.LeastSquares import leastSquaresFit as LS
import pyfits
import pylab

__all__=['run_clump','Av_d','convert']
print """Red clump method to be supreseeded by one base on 
histogramdd"""

histmin=0.5
histmax=2.5
binsize=0.1
photcut="B"

def model(p,x):
    """returns a power-law plus Gaussian feature"""
    norm,index,x0,sig,peak = p
    return abs(norm) * x**index + abs(peak)*exp(-((x-x0)**2)/(2*sig**2))

def fit_hist(bin,count,normguess=None,indexguess=-1,x0guess=1,sigguess=0.2,
             peakguess=None):
    """fit a histogram to the model above"""
    data = zeros((len(bin),3))*1.0
    data[:,0] = bin
    data[:,1] = count
    data[:,2] = 1 # pylab.clip(sqrt(count),1,1e8)
    normguess = normguess or average(count)
    peakguess = peakguess or normguess/10
    return LS(model, (normguess,indexguess,x0guess,sigguess,peakguess), data)

def run_clump(J,K,flags,startslices=arange(8,14,.2),endslices=arange(9,15,.2)):
    """ fit 2mass data J and K magnitudes to give Av-d graph by fitting the red clump
    J and K are apparent magnitudes from 2mass, flags is 2mass photometry flags "Qflg"
    Fit returns set of J-K of clump, Kvalues, sigmas, and number of stars in each clump"""
    J_Klist=[]  # output lists
    myhistmin=histmin
    Klist = []
    sigs=[]
    nstars=[]
    p4=[None,1e10]
    p1=[None,1e10]
    p2=[None,1e10]
    p3=[None,1e10]
    if len(startslices) != len(endslices):
        print "Must have same number of start and end points!"
        return
    global bins,counts,histmax
    for index,flag in enumerate(flags):  #reject bad photometry
        if (flag[0]>photcut)or(flag[2]>photcut):
            J[index]=0
    K = K[pylab.where(J>0)]
    J = J[pylab.where(J>0)]
    pylab.clf()
    pylab.plot(J-K,K,'k.')
    pylab.axis([0,2.5,15,5])
    pylab.xlabel('$J-K$')
    pylab.ylabel('$K$')
    print len(J), "valid stars for analysis"
    col = J-K
    x0guess=0.75
    sigguess=0.1
    norm = 1
    indie = 0
    peakie=4
    pylab.figure(2)
    iterator = range(len(startslices))
    for i in iterator:
        this_slice = endslices[i]
        prev_slice = startslices[i]
        colext = col[pylab.where((K>prev_slice)*(K<this_slice))]
        print "Slice %f -> %f  has %i stars"%(prev_slice,this_slice,len(colext))
        pylab.clf()  
        colext2=colext[pylab.where((colext>0.5)*(colext<histmax))]
        counts,bins,patches = pylab.hist(colext2,bins=arange(myhistmin,histmax,binsize),facecolor=(.5,.5,1))
        bins += binsize/2 #middle of hist bins
        fail=0
        try: p1 = fit_hist(bins,counts,x0guess=x0guess,sigguess=sigguess,normguess=norm,
                     indexguess=indie,peakguess=peakie)
        except: fail+=1 #try other guesses too
        try: p2 = fit_hist(bins,counts,x0guess=x0guess+0.2,sigguess=sigguess,normguess=norm,
                     indexguess=indie,peakguess=peakie)
        except: fail+=1
        try: p3 = fit_hist(bins,counts,x0guess=x0guess+0.4,sigguess=sigguess,normguess=norm,
                     indexguess=indie,peakguess=peakie)
        except: fail+=1
        try: p4 = fit_hist(bins,counts,x0guess=x0guess+0.6,sigguess=sigguess,normguess=norm,
                     indexguess=indie,peakguess=peakie)
        except: fail+=1
        try: p5 = fit_hist(bins,counts,x0guess=x0guess-0.2,sigguess=sigguess,normguess=norm,
                     indexguess=indie,peakguess=peakie)
        except: fail+=1
        bestfit = pylab.array([p1[1],p2[1],p3[1],p4[1]]).argmin()
        p = (p1,p2,p3,p4)[bestfit] #pick the fit with best chi^2
        p1=[None,1e10]
        p2=[None,1e10]
        p3=[None,1e10]
        p4=[None,1e10]
        p5=[None,1e10]
        if fail==5: #none of the fits worked
            iterator.insert(0,0) #cheat to repeat this phase
            x0guess = input("Fitting failed - enter new x0guess: ")
            sigguess= input("New sig: ")
            indie = input("New index: ")
            peakie = input("New peak: ")
            norm = input("New norm: ")
            continue
        norm,indie,x0guess,sigguess,peakie = p[0]
        try:
            pylab.plot(bins, abs(p[0][0]) * bins**p[0][1] + exp(-(bins-p[0][2])**2/p[0][3]**2
                                                 ) * abs(p[0][4]), 'go')
            bins2=arange(min(bins),max(bins),0.01)
            pylab.plot(bins2, abs(p[0][0]) * bins2**p[0][1] + exp(-(bins2-p[0][2])**2/p[0][3]**2
                                                 ) * abs(p[0][4]), 'g-')
        except:
            pass
        s=raw_input("<r> to repeat, <s> to shift, <q> to quit: ")
        if s=='r':
            iterator.insert(0,0) #cheat to repeat this phase
            x0guess = input("Repeating: Enter new x0guess: ")
            sigguess= input("New sig: ")
            indie = input("New index: ")
            peakie = input("New peak: ")
            norm = input("New norm: ")
            continue
        if s=='q':
            break
        if s=='s':
            myhistmin+=0.2
        Klist.append((this_slice+prev_slice)/2)
        J_Klist.append(x0guess)
        sigs.append(sigguess)
        nstars.append(17.7*sigguess*peakie)
        prev_slice=this_slice
        sigguess = min((max((sigguess,0.15)),2.45))
        x0guess = min((max((x0guess,0.8)),2))
    return pylab.array(J_Klist),pylab.array(Klist),pylab.array(sigs),\
           abs(pylab.array(nstars))

def Av_d(myfile,maxstar=None):
    """run red-clump fitting on the photometry in file, being a 2mass table in a particular
    direction, either FITS or ascii"""
    if myfile[-4:]=="fits" or myfile[-3:]=='fit' or myfile[-3:]=="FTS":
        fred = pyfits.open(myfile)
        pylab.figure(1)
        pylab.clf()
        table = fred[1].data
        fred.close()
        J = table.field("Jmag")
        K = table.field("Kmag")
        flags = table.field("Qflg")
        r = table.field("_r")
        if maxstar:
            print "Max radius: %f arcmin"%r[maxstar]
        else: print "Max radius: %f arcmin"%max(r)
    else:   #assume ascii file, with columns id,ra,dec,j,jerr,h,herr,k,kerr,flag...
        fred = open(myfile)
        mylines = fred.readlines()
        J=[]
        K=[]
        flags=[]
        for aline in mylines:
            words = aline.split()
            if len(words)<10: continue
            J.append(float(words[3]))
            K.append(float(words[7]))
            flags.append(words[9])
        J=pylab.array(J)
        K=pylab.array(K)
    if maxstar:
        J = J[:maxstar]
        K = K[:maxstar]
        flags=flags[:maxstar]
    fred.close()
    pylab.figure(1)
    pylab.clf()
    pylab.plot(J-K,K,'k.')
    pylab.xlabel("$J-K$")
    pylab.ylabel("$K$")
    first = eval(raw_input('Slice star values (python syntax): '))
    last = eval(raw_input('Slice end values (python syntax): '))
    jk,k,sig,nstar = run_clump(J,K,flags,first,last)
    pylab.figure(1)
    sam = pylab.errorbar(jk,k,yerr=(first-last)/2,xerr=sig/pylab.sqrt(nstar),capsize=0)
    pylab.setp(sam[0],marker='d',mec='w',mfc='w',markersize=5,ls='None')
    pylab.setp(sam[1],c='w',lw=2)
    sam = pylab.errorbar(jk,k,yerr=(first-last)/2,xerr=sig/pylab.sqrt(nstar),capsize=0)
    pylab.setp(sam[0],marker='d',mec='c',mfc='c',markersize=4,ls='None')
    pylab.setp(sam[1],c='c',lw=1)
    pylab.axis([0,2.5,15,5])
    av,d,err = convert(jk,k,sig,nstar)
    pylab.figure(2)
    pylab.clf()
    pylab.errorbar(d,av,err)
    return jk,k,sig,nstar
    
def convert(jk,k,sig,nstar):
    """take the outputs of the above, and return A_V and d with
    simple errors assuming good gaussian errors, using numbers
    from Drimmel et al. (2003)amd Schlegel et al. (1998)"""
    Av = (jk - 0.75)/0.164
    d = 10**(0.2*(k+1.65-0.112*Av))*10
    # note correction for reddening
    err = sig/(pylab.sqrt(nstar)*0.164)
    return Av,d,err
