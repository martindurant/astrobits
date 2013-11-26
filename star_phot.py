#!/usr/bin/python
import sys
try:
    from numpy import r_,array,log10,zeros,where,unravel_index,nan,load,newaxis,sqrt,std,isfinite
except:
    print """
Requires numpy to run. This can be installed individually
(e.g. "sudo apt-get install python-numpy" on a debian system)
or with a package like scisoft, enthought, python(x,y)...
"""
    raise

__doc__ = """
STAR CLASSIFIER
A small program to give you a rough spectral class based on broad-band photometry.
Input UBVRIJHK magnitudes are compared to a library of models computed using
pysynphot from the k93models (Kurutz, 1993, see
http://www.stsci.edu/hst/HST_overview/documents/synphot/AppA_Catalogs9.html#49686).
It will return the effective temperature, metallicity, gravity and reddening of the
nearest matching model.

Usage:
star_phot.py U B V R I J H K
or
star_phot.py U eU B eB V eV R eR I eI J eJ K eK

where the capital letter parameters are numerical magnitudes in those standard
bands, and "e" parameters are uncertainties. To exclude any particular band from
the fit, give its magnitude as "nan" or make the uncertainty negative.

The provided grid of magnitudes gives a reasonable trade-off between size and
coarseness. If you have pysynphot installed, you can make your own with
make_grid(), choosing the parameters and model catalogue to use.

I hope you find this useful, but remember, it is only meant as a quick-look tool,
not for serious categorisation.
Martin Durant
"""

Teff = r_[3000:10000:250,10000:13000:500,13000:35000:1000,35000:50000:2500]
logZ = r_[-5:-0.49:0.5,-0.3:0.31:0.1,0.5,1]
logg = r_[0:5.1:0.5]
Eb_v = r_[0:2:0.05]

def load_grid(filename):
    """Input grid and values from a file previously made with make_grid and saved with numpy.savez"""
    fred = load(filename)
    return fred['grid'],fred['Teff'],fred['logZ'],fred['logg'],fred['Eb_v']

def find_mags(T,logZ,logg,E,cat='k93models'):
    """Consult k93models for the given parameters, and return absolute mags.
WARNING: values for T<3500K are ill-defined."""
    import pysynphot as S
    mags = [S.ObsBandpass(i) for i in ['U','B',"V","R","I","J","H","K"]]
    model = S.Icat(cat,T,logZ,logg)
    model2= model*S.Extinction(E)
    out = zeros(8.)
    for i in range(8):
        try:
            out[i] = S.Observation(model2,mags[i]).effstim('vegamag')
        except:
            out[i] = nan
    return out

def weighted_mean(x,dx):
    """Mean of x weighted by the uncertainties in dx. Returns mean
and uncertainty (assumes Gaussian uncertainties)."""
    shape = x.shape
    ind = isfinite(x).ravel()
    num = (x.ravel()[ind]/dx.ravel()[ind]**2)
    num2=num.reshape(shape)
    den = (1/dx.ravel()[ind]**2)
    den2=den.reshape(shape)
    mean = num2.sum(axis=4)/den2.sum(axis=4)
    return mean

def find_star(grid,mags,err=1,Teff=Teff,logZ=logZ,logg=logg,Eb_v=Eb_v):
    """Given a grid calculated by make_grid, find the best-matching set of magnitudes to the input magnitudes.
If you only have a subset of UBVRIJHK, either set those errors massive,
or input magnitudes as nan."""
    diff = grid - mags
    shape=diff.shape
    if (type(err) is int) or (type(err) is float):
        dx = diff-diff+err #if err is just a number
    else:
        dx = diff-diff
        dx[...] = err[newaxis,newaxis,newaxis,newaxis,:]
    dx2 = dx.copy()
    dx[grid>99] *= 1e6 # Do not use non-values in the grid for finding the mean offset
    mean = weighted_mean(diff[...,isfinite(mags)],dx[...,isfinite(mags)]) #mean offset
    ind = isfinite(mags)
    err = std(((diff.transpose()-mean.transpose())/dx2.transpose()).transpose()[:,:,:,:,ind],axis=4)
    err[grid.min(axis=4)>99] = 1e10
    ind = err.argmin()
    ind = unravel_index(ind,err.shape)
    print "Catalogue index:" ,ind
    print "Chi2: ",err.min(), "  Magnitude offset: ",-mean[ind]
    print "Delta-mags:", diff[ind]-mean[ind]
    return Teff[ind[0]], logZ[ind[1]], logg[ind[2]], Eb_v[ind[3]]

def make_grid(Teff=Teff,logZ=logZ,logg=logg,Eb_v=Eb_v,cat='k93models'):
    """ For each combination of T, logZ, logg, E(B-V), produce
absolute BVRIJHK mags. Takes a while."""
    import pysynphot as S
    U = S.ObsBandpass('U')
    B = S.ObsBandpass('B')
    V = S.ObsBandpass('V')
    R = S.ObsBandpass('R')
    I = S.ObsBandpass('I')
    J = S.ObsBandpass('J')
    H = S.ObsBandpass('H')
    K = S.ObsBandpass('K')
    mags = [U,B,V,R,I,J,H,K]
    effs = array([mag.avgwave() for mag in mags])
    fred = open('out.text','w')
    outarray = zeros((len(Teff),len(logZ),len(logg),len(Eb_v),8))*1.
    for i in range(len(Teff)):
        fred.write("%e\n"%Teff[i])
        for j in range(len(logZ)):
            for k in range(len(logg)):
                for l in range(len(Eb_v)):
                    outarray[i,j,k,l,:] = find_mags(Teff[i],logZ[j],logg[k],Eb_v[l],cat=cat)
    return where(outarray<0,outarray,99.9)

mygrid = 'grid.npz'

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv)<2:
        print __doc__
        return 0
    try:
        inmags = array(argv[1:]).astype(float)
    except:
        print "Expect numerical entry only"
        return 2
    if len(inmags)==8:
        mymags = inmags
        err=1
    elif len(inmags)==16:
        mymags = inmags[::2]
        err = inmags[1::2]
        mymags[err<0] = nan
    else:
        print """Enter either eight values (UBVRIJHK magnitudes) or sixteen values
(Umag Uerr ...). Where neasurement does not exist, enter "nan" for the
magnitude or a negative number for the uncertainty."""
        return 2   
    grid,Teff,logZ,logg,Eb_v = load_grid(mygrid)
    Teff,logZ,logg,Eb_v=find_star(grid,mymags,err,Teff,logZ,logg,Eb_v)
    print "Teff: ", Teff, "  log Z: ",logZ
    print "log g:", logg, "  E(B-V):",Eb_v
    return 0
    
   

if __name__=='__main__':
    print "Star Classifier, Martin Durant 2010\n"
    sys.exit(main())
