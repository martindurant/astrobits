from pyraf import iraf
iraf.martin()
import pyfits
from numpy import mean,median
import os
from glob import glob

def irstack(myfiles):
    data = []
    tempfiles = glob("templist")+glob("*b.fits")+glob("*r.fits")+glob("flat.fits")
    for myfile in tempfiles:
        os.remove(myfile)
    for myfile in myfiles:
        data.append(pyfits.getdata(myfile))
    #sky subtract and replace with median
    mediansky = mean([median(data[i]) for i in range(len(data))])
    for i in range(len(myfiles)):
        fred = pyfits.open(myfiles[i])
        im = fred[0].data
        im2 = (im.transpose() - median(im,axis=0)).transpose() + mediansky
        fred[0].data = im2
        fred.writeto("%ibb.fits"%i,'ignore',True)
    iraf.imcomb("*bb.fits","flat",combine="median",reject="sigclip",lsigma=3,hsigma=2)
    flat = pyfits.getdata('flat.fits')
    flat /= median(flat)
    fred.writeto('flat.fits','ignore',True)
    for i in range(len(myfiles)):
        fred = pyfits.open(myfiles[i])
        im = fred[0].data
        im2 = ((im/flat).transpose() - median(im,axis=0)).transpose()
        fred[0].data = im2
        fred.writeto("%irb.fits"%i,'ignore',True)
    iraf.files("*rb.fits",Stdout="templist")
    iraf.stack("templist","output",1,'none')
