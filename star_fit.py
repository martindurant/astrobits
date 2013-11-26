from pylab import *
from Numeric import *
import scipy.optimize
from scipy.optimize import leastsq as leastsq

def gauss(amp,x0,sigma,x):
    # gaussian which works on either scalar or array
    return amp * exp( - ((x-x0) / sigma)**2)

def gauss2darr(zp,amp,x0,y0,sigmax,sigmay,x,y):
    # for 2d gaussian with two 1d array inputs
    # amp is total amplitude
    return  gauss(amp,x0,sigmax,x) * gauss(1,y0,sigmay,y)[:,NewAxis] + zp

def gauss2d(amp,x0,y0,sigmax,sigmay,x,y):
    return amp * exp( - ((x-x0) / sigmax)**2 ) * exp( - ((y-y0) / sigmay)**2 )

def residuals(p,z,y,x):
    zp,amp,x0,y0,sigmax,sigmay = p
    out = (z-gauss2darr(zp,amp,x0,y0,sigmax,sigmay,x,y))
    out.shape=(len(x)*len(y),)
    return out

def relative(p,z,x0,y0,sigmax,sigmay,y,x):
    amp,zp=p
    # here, only 'amp' and 'zp' is allowed to vary
    out = (z-gauss2darr(zp,amp,x0,y0,sigmax,sigmay,x,y))
    out.shape=(len(x)*len(y),)
    return out

def relative_2star(p,z,xA,yA,xB,yB,sigmax,sigmay,y,x):
    ampA,ampB,zp=p
     # here, only 'ampA', 'ampB' and 'zp' is allowed to vary
    out = (z-gauss2darr(zp,ampA,xA,yA,sigmax,sigmay,x,y)\
           -gauss2darr(zp,ampB,xB,yB,sigmax,sigmay,x,y))
    out.shape=(len(x)*len(y),)
    return out

def fullfit(image,zpi,ampi,xi,yi,sigxi,sigyi):
    # fit for general Gaussian
    x,y = image.shape
    x = arange(x)
    y = arange(y)
    fred = leastsq(residuals,(zpi,ampi,xi,yi,sigxi,sigyi),args=(image,y,x))
    return fred[0]  # contains zp,amp,x,y,sigx,sigy

def relfit(image,zpi,ampi,xi,yi,sigxi,sigyi):
    # fit for amplitude only
    x,y = image.shape
    x = arange(x)
    y = arange(y)
    fred = leastsq(relative,(ampi,zpi),args=(image,xi,yi,sigxi,sigyi,y,x))
    return fred[0]

def relfit_2star(image,zpi,ampA,ampB,xA,yA,xB,yB,sigxi,sigyi):
    # fit for amplitude only
    x,y = image.shape
    x = arange(x)
    y = arange(y)
    fred = leastsq(relative_2star,(ampA,ampB,zpi),args=(image,xA,yA,xB,yB,sigxi,sigyi,y,x))
    return fred[0]

