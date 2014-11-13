from pylab import median,diff,std,sqrt

def EW(wave,flux,a,b,c,d,e,f):
    """calculate the equivalent width of a feature in a spectrum.
    The letters deffine three regions in wavelength, two for the
    continuum and one for the feature. The output EW is in the units
    of wave.
    Returns EW, error, continuum and std"""
    ind = (wave>a)*(wave<b) + (wave>e)*(wave<f)
    sig = std(flux[ind])
    cont = median(flux[ind])
    ind = (wave>c)*(wave<d)
    deff = (sum(flux[ind]) - sum(ind)*cont)/cont
    err = sig*sqrt(sum(ind))/cont
    ew = deff*median(diff(wave[ind]))
    return ew,abs(ew*err),cont,sig

#to do: graphical version of this? Full spectral suite???
