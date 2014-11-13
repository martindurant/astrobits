"""Fits real periodogram data, assuming that the powers
follow a chi2 distributioni around a model (power-law plus
constant, typically). This is a real statistical analysis,
since least-squares fitting assumes gaussian statistics. The
chi2(2) distribution is strictly valid only for Fourier
periodograms."""

import scipy

from scipy.stats import chi2
from scipy.optimize import fmin,leastsq
from numpy import multiply,arange,log


def fitfunc(p,x):
    """model to fit for, parameters in p. Here a power-law
plus constant, but could be anything else"""
    return p[1]*x**p[2] + p[0]

def errfunc(p,x,y):
    """produce goodness measures for fitfunc applied to independant
variable x compared with data y for chi2(2) statistics"""
    chis = 2*y / fitfunc(p,x) #values relative to current model
    probs = chi2.pdf(chis,2) #probability for each value
    return log(probs).sum() #function goodness

def fullfit(freq,power,f0=0):
    """Perform chi2-based fitting of power spectrum.
f: frequencies
p: corresponding powers
f0: lowest frequency to use"""
    assert len(freq)==len(power)
    fnew = freq[freq>f0]
    pnew = power[freq>f0]
    p0 = [pnew.mean(),pnew[0]*fnew[0],-1]
    return fmin(errfunc,p0[:],args=(fnew,pnew))
