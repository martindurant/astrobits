"""Misc astronomical data"""
from startup import *

#photometric
filters = {"K":0, "H":1, "J":2, "I":3, "R":4, "V":5, "B":6}
red = array([0.112,0.176,0.276,0.594,0.819,1.015,1.321]) #A/A_V
l = array([2.19,1.63,1.22,0.797,0.638,.545,.436])*1e-6 #central wavelength (m)
zp = array([.636,1.02,1.57,2.42,3.06,3.60,4.00])*1e-20

redfunc = lambda nu: 1.07186134e-01-1.14018516e-15*nu+1.00929467e-29*nu**2-1.13857192e-44*nu**3+4.15062871e-60*nu**4
#polynomeal fit to reddening curve

c=299792458. #m/s
h=6.62606896e-34 #Js
k=1.3806504e-23 #J/K

def X_ext(E, nH,dE=0):
    """Extinction from XSPEC 'wabs', for input energy E
in keV and N_H (10^24 cm-2)returns factor (<1) """
    energy, absorb, another = loadtxt('/usr/local/heasoft-6.11/spectral/modelData/absorp.dat',unpack=True)
    if dE is 0:tau = lininterp(energy,absorb,E)
    else:
        tau = E-E
        for i in range(len(E)):
            tau[i] = lininterp(energy,absorb,linspace(E[i]-dE[i],E[i]+dE[i],100)).mean()
    correct = exp(-nH*tau)
    return correct

def bb(T,x,lam=True,R=1,D=1):
    """Black body flux F_lam(lam) or F_nu(nu); calibrate for Radius R km
and distance D pc. For x in m or Hertz, gives back erg."""
    if lam:
        F=2*h*c**2/x**5/(exp((h*c)/(x*k*T))-1)
    else:
        F=2*h*x**3/c**2/(exp((h*x)/(k*T))-1)
    return (R*1e3)**2/(D*30.857e15)**2*F*1e7

def red_nu(nu,ebv):
    """Calls reddening below for nu (in Hz) and returns an extinction (number<1)
factor"""
    lam = c/nu
    red = reddening(lam*1e6)
    return 10**(-0.4*red*ebv)

def reddening(lam):
    """Seaton (1979, MNRAS, 187, 73)'s reddening function deep into the UV.
Input lambda in um, output A_lam/E(B-V). Returns zero where not defined.
Assumes R=Av/E(B-V)=3.2 ."""
    y = 1/lam
    out = where((y<=10)*(y>=7.14),16.17-3.2*y+0.2975*y**2,0)
    out = where((y<=7.14)*(y>=3.65),2.29+0.848*y+1.01/((y-4.60)**2+0.280)  ,out)
    out = where((y<=3.65)*(y>=2.7),1.56+1.048*y+1.01/((y-4.60)**2+0.280)  ,out)
    out = where((y<=2.7)*(y>=1.1),3.1*(1 + 0.17699 * (y-1.82) - 0.50447 * (y-1.82)**2.0 -0.02427 * (y-1.82)**3.0 + 0.72085 * (y-1.82)**4.0 \
           +0.01979 * (y-1.82)**5.0 - 0.77530 * (y-1.82)**6.0 +0.32999 * (y-1.82)**7.0) + (1.41338 * (y-1.82) +2.28305 * (y-1.82)**2.0 \
           +1.07233 * (y-1.82)**3.0 - 5.38434 * (y-1.82)**4.0 \
           -0.62251 * (y-1.82)**5.0 + 5.30260 * (y-1.82)**6.0 -2.09002 * (y-1.82)**7.0)     ,out)
    out = where((y<=1.1)*(y>=0.3),3.1* (0.574 * y**1.61) + (-0.527 * y**1.61),  out)
    return out

#phot zero point (erg/s/cm^2/Hz)
def flux(mag,band):
    """get f_nu flux and freq for mag in filter 'band'.
    Flux is in erg/s/cm2, nu in Hz"""
    zero = zp[filters[band]]
    flux = zero * 10**(-0.4*mag)
    nu = c/l[filters[band]]
    return flux,nu

def calib_flux(target,standard,table,scale=1,interp='lin',window=5):
    """For target and standards extracted spectra and table tabulated values,
    return (wav,flux) calibrated spectrum. This is done by linear interpolation
    only (use interp!='lin' for splines).
    Scale is incase the fluxes are tabulated in odd units.
    Window is for smoothing the response function, with a window of this many wav bins.
    Also returns the response function for checking."""
    wav,pix = wavelength_load(standard)
    wav2,pix2=wavelength_load(target)
    lam,flux=load(table,usecols=[0,1],unpack=True)
    resp = make_resp(lam,flux,wav,pix)
    if interp!='lin':
        resp2=splineinterp(lam,resp,wav2)
    else:
        resp2 = lininterp(lam,resp,wav2)
    resp2 = where(resp2>0,resp2,0)
    if smooth>1:
        resp2 = smooth(resp2,window)
    scale = diff(lam).mean()/diff(wav).mean()
    print "Wavelength rescale: ",scale
    return wav2,pix2*resp2/scale,resp2/scale
    

def make_resp(lam,flux,wav,pix,binsize=None):
    """(lam,flux) being tabulated standard data, find the average response
    at each lam given measured (wav,pix). You probably then want to interpolate
    this to your target spectrum using myinterp. Response function is for
    multiplying into count values."""
    resp = flux-flux
    binsize = binsize or diff(lam).mean()
    for i in range(len(lam)):
        ind = (wav>=lam[i]-binsize/2.)*(wav<lam[i]+binsize/2.)
        resp[i] = flux[i]/pix[ind].mean()
    return resp

def wavelength_load(fitsfile):
    """For an IRAF-generated spectrum with wavelength scale information in the
    fits header, return wavelength versus counts"""
    fred = pyfits.open(fitsfile)
    pix = fred[0].data
    if len(pix.shape)>1: #"full output"
        pix = pix[0,0]
    try:
        wav = arange(len(pix))*fred[0].header['CDELT1'] + fred[0].header['CRVAL1']
    except:
        wav = arange(len(pix))*fred[0].header['CD1_1'] + fred[0].header['CRVAL1']
    return wav,pix

c=2.9979e8   #mks
h=6.6261e-34 #mks

mag0142 = array([ 20.0,  21.1 ,  22.3,  23.84,  24.89,  25.62,  28.1 ])
err0142 = array([ 0.1,  0.1,  0.1 ,  0.04,  0.05,  0.08,  0.3 ])
# Israel et al (2004) for K'HJ, Hulleman et al (2004) for IRVB
                
#mag0142 = array([ 19.85,  20.6 ,  22.07,  23.84,  24.89,  25.62,  28.1 ])
#err0142 = array([ 0.11,  0.08,  0.1 ,  0.04,  0.05,  0.08,  0.3 ])
# Morii et al (2005) for K'HJ, Hulleman et al (2004) for IRVB

def rotationcurve(l,b,d=arange(0.1,25,0.1)):
    """Calculate rotation curve through the galactic plane in
    the direction longitude, l, and latitude, b (galactic, degrees).
    Returns d in kpc and v in km/s"""
    theta0 = 220  #local circular velocity
    r0 = 8.5      #distance to the Galactic centre
    dproj = d*cos(b/180.*3.14159)
    r=sqrt(r0**2+dproj**2-2*r0*dproj*cos(l/180.*3.14159))
    w=1.00746*theta0/r-0.017112*theta0/r0
    v=(w-theta0/r0)*r0*sin(l/180.*3.14159)*cos(b/180.*3.14159)
    return d,v

def hms(deg,astext=1):
    h=int(deg)
    mfrac = (deg-h)*60*sign(deg)
    m=int(mfrac)
    s = (mfrac-m)*60
    if astext:
        return "%i:%02i:%5.3f"%(h,m,s)
    return h,m,s
