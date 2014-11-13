import pyfits
from numpy import histogram2d,argmax,array,arange,floor,inf
import ndimage

def image(infile,outfile,xcol,ycol,smoothing=None,ext=0):
    """Read FITS event data, and produce FITS image by histogram, copying the WCS
from the original header. xcol and ycol are the variables to use, either
column index or name. Filter before applying this!"""
    myfile = pyfits.open(infile)
    head = myfile[ext].header
    data = myfile[ext].data
    myfile.close()
    if xcol.__class__ is " ".__class__  :
        xcol = argmax(array(data.names) == xcol)
    if ycol.__class__ is " ".__class__  :
        ycol = argmax(array(data.names) == ycol)
    x = data.field(xcol)
    xmin = floor(x.min())
    x -=x.min()
    y = data.field(ycol)
    ymin = floor(y.min())
    y -=y.min()
    cols = (xcol,ycol)
    mins = (xmin,ymin)
    im,Y,X = histogram2d(y,x,bins=(arange(y.max()+1),arange(x.max()+1)))
    if smoothing:
        im = ndimage.gaussian_filter(im,smoothing,mode='nearest')
    hdu = pyfits.ImageHDU(im)
    wordsin = ['TCTYP','TCRVL','TCRPX','TCDLT','TCUNI']
    wordsout= ['CTYPE','CRVAL','CRPIX','CDELT','CUNIT']
    for i in range(1,3):
        for j in range(len(wordsin)):
            hdu.header.update(wordsout[j]+str(i),head[wordsin[j]+str(cols[i-1]+1)])
        hdu.header.update('CRPIX'+str(i),head['TCRPX'+str(cols[i-1]+1)]-mins[i-1])
    hdu.writeto(outfile)
