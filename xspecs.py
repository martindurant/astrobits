import pyfits
from pylab import *
from startup import weighted_mean as WM

def data_in(file):
    """read XMM spectra file with columns CHANNEL,FLUX,ERROR"""
    fred = pyfits.open(file)
    table = fred[1].data
    fred.close()
    wave = table.field("CHANNEL")
    flux_l = table.field("FLUX") 
    err = table.field("ERROR")
    err = where(((flux_l<=0)+(flux_l>0)),err,zeros(len(flux_l))+10) 
    flux_l = where(((flux_l<=0)+(flux_l>0)),flux_l,zeros(len(flux_l)))
    return array(wave),array(flux_l),array(err)

def constant_rebin(wavein,fluxin,errin,factor):
    """rebin by a constant factor, using wighted means"""
    if (len(wavein)<>len(fluxin))&(len(wavein)<>len(errin)):
        print "\nArrays not of same size!"
        return
    binsin = len(wavein)
    binsout = int(binsin/factor)
    waveout = arange(binsout)*0.1
    fluxout = arange(binsout)*0.1
    errout  = arange(binsout)*0.1
    for i in range(binsout):
        startbin = i*factor
        endbin = startbin+factor
        waveout[i] = mean(wavein[startbin:endbin])
        fluxout[i],errout[i] = WM(fluxin[startbin:endbin],errin[startbin:endbin])
    return waveout,fluxout,errout

def adaptive_rebin(wavein,fluxin,errin,relerr,minerr=0,minrelerr=0):
    """
    rebin the fluxes into wavebins with errors below relerr.
    less than minerr (absolute) or minrelerr (relative) always
    gets a bin to itself.
    """
    if (len(wavein)<>len(fluxin))&(len(wavein)<>len(errin)):
        print "arrays not of same size"
        return
    binsin = len(wavein)
    waveout=[]
    fluxout=[]
    errout=[]
    index=0
    while index < binsin-1:
        start = index
        end = index+1
        try:
            err = errin[index]
            total = fluxin[index]
            myrelerr = err/total
        except:
            myrelerr=1e9
        while (abs(myrelerr) > relerr):
            end=end+1
            N=end-start
            total = sum(fluxin[start:end])/N
            err = 1/sqrt(sum(1/errin[start:end]**2))
            try:
                myrelerr = err/total
            except:
                myrelerr = 1e9
            if end == binsin:
                end=binsin-1
                break
            if errin[end]<minerr or abs(errin[end]/fluxin[end])<minrelerr:
                if N==1:
                    break
                else:
                    end = end-1
                    break
        fluxout.append(total)
        errout.append(err)
        waveout.append((wavein[end-1]+wavein[start])/2)
        index=end
    return array(waveout),array(fluxout),array(errout)

k=1.38066e-23 #Botzmann
h=6.62608e-34 #Planck
c=2.99792e8 #light

def blackbody(l,T):
    """for the wavelengths l (m), return the blackbody
    curve deffined by temperature T (in Kelvin).
    Returns units J/s/m2/m/sr (i.e., energy)"""
    I = (2*h*c**2/l**5) / (exp(h*c/(l*k*T))-1) 
    return I 

def bb(e,kT):
    """as blackbody, but takes e in keV, kT in keV and returns
    units erg/s/cm2/Hz/sr"""
    nu = (e * 1.602e-16)/6.626e-34
    I =  (2*h*(nu**3)/c**2) / (exp(e/kT)-1)
    # in J/s/m2/Hz/sr
    return I * 1e7*1e-4  

def luminosity(r,kT):
    """bolometric luminosity for a blackbody kT(keV) at d(pc) and
    radius r(km) in erg/s.
    SHOULD'T THIS JUST BE STEFAN-BOLTZMANN??"""
    e=arange(.001,100,0.001)
    I=bb(e,kT)
    return sum(I * 4*pi*(r*1e5)**2 * 2)*0.001 * 0.92e17
