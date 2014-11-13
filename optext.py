__doc__="""Optimal photometry, for one faint star and one bright reference star.
Optimised for timeseries creation. Data variables are module wide.
optext: runs the routine"""

from Scientific.Functions.LeastSquares import leastSquaresFit as LS


#fitting data (used but not changed in functions)
FWHM = 3.
inannmult = 2. #normaly ~twice FWHM
inannguess = 6.
outannmult = 4.5. #should have few hundred sky pixels
outannguess = 13.
#NB: make sure if using a median-flat, that it is made of more images
#than the number of sky pixels
skyguess = 10000.
skyvar = 200.
relx = -20.
rely = 35.
sigmas = 5

def 2dgauss(height,x,y,width1,width2,posangle,inannulus):
    """returns a model psf, whose size is no more than inannulus"""
    return
    
def skewgauss(height,x,width,skew):
    """returns a skewed gaussian for fitting the sky - see Irwin(1995)"""
    return

def fit_reference(image, sky, skyvar, xguess, yguess):
    """Fit 2D elliptical gaussian to image, using *guess as initial
    parameters and ignoring all values sigmas number of skyvar above
    the sky or less.
    Returns position x,y, perpendicular widths width1,width2 and
    position angle"""
    return

def fit_sky(image, x, y, inannulus, outannulus):
    """Fit skewed Gaussian distribution to the sky values within
    the annuli centred on x,y in image.
    Returns sky value and varience."""
    return

def fit_object(image, x, y, width1, width2, posangle, sky):
    """Find the flux in the object at x,y using the gaussian shape deffined.
    Returns flux and error"""
    return

def optext(image, refx, refy):
    """Use Tim Naylor's formulism for optimal photometry. First, fits sky
    around bright reference star, then finds the star shape, then fits for
    the faint star of interest. Uses 2 iterations for the bright star, first
    to get rough parameters which set the annuli for the sky, then with the
    imporoved sky value."""
    brightsky,brightvar = fit_sky(image, refx, refy, inannguess, outannguess)
    #fit rough sky with typical values
    x,y,width1,width2,posangle = fit_reference(image,brightsky,brightskyvar,refx,refy)
    # get rough fit
    inannulus = max((width1,width2)) * inannmult
    outannulus = max((width1,width2)) * outannmult
    brightsky,brightvar = fit_sky(image, x, y, inannulus, outannulus)
    # refit sky around reference
    x,y,width1,width2,posangle = fit_reference(image,brightsky,brightskyvar,x,y)
    # good fit for reference
    faintsky,faintvar = fit_sky(image, x+relx, y+rely, inannulus, outannulus)
    # fit sky for faint object
    flux, error = fit_object(image, x+relx, y+rely, width1, width2, posangle, faintsky)
    # fit faint object
    return flux, error
