from startup import prctile,zeros,diff,arange,std,rebin,median
from scipy.stats import distributions
rvs = distributions.norm.rvs

def do(func,p,err,E,prct=[16,50,84],num=1000,**kwargs):
    """Take a function, and calculate the mean flux and +/-1sig
(or prct) values of the flux in the given energy bins.
Works by monte-carlo random sampling of the parameters,
assuming they are Gaussian.
[Assumes Ebins are regular, and 10 flux bins per output bin]
func(p,x,*args,**kwargs) must give a flux for parameters p,
and energy/wave x"""
    flux = zeros((num,len(E)))
    for i in range(num):
        p1 = [rvs(p[j],err[j]) for j in range(len(p))]
        out = func(p1,E,**kwargs)
        flux[i] = out
    centre = flux.mean(axis=0)
    lower = centre-centre
    upper = centre-centre
    for i in range(len(E)):
        lower[i],centre[i],upper[i] = prctile(flux[:,i],prct)
    return lower,centre,upper
