"""Take a list of files and known star coordinates, and
perform photometry on them all, either with apertures (phot)
or by PSF fitting (daophot, which required additional
parameters and is apropriate to poor S/N or crowded fields).
Makes extensive use of iraf tasks; set all photometry parameters
before running:
datapars - for data characteristics
centerpars - finding the reference star on each image.
centerpars, photpars, fitskypars - for controling aperture photometry
daopars - for controling daophot

filelist: set of image files, in IRAF syntax (image.fits[1][*,*,2] etc);
can be more than one per cube.
coords: name a file containing all star coords for photometry, based on
an image unshifted relative to (0,0) in the shifts list. Be pure numbers
for phot method, .mag or .als for daophot method.
shifts: name a file containing shifts, a tuple of shifts arrays, image
header keywords (tuple of two= or None for no shifts
refstar: coords of star for deriving (x,y) offset, as in coords
timestamp: source of the timing information: a header keyword, delta-t
for uniform sampling or a file with times (in whatever formate you'll be
using later.

psf: whether to use daophot or aperture phot for analysis. If this is a
filename, that is the PSF profile to use for every image; if it is "True",
make a new PSF for every image. Pars below only for full PSF fitting
pststars: a .pst file from daophot, listing the IDs of stars for making
the PSF for each image. NB: DAOphot refuses to measure any star with SNR<2.

ids: which stars are interesting, by ID (in input coord list order)
coords: starting well-measured coords (pdump-ed from a .als, perhaps).
"""
import os
import numpy
from glob import glob
import pyfits
from pylab import find
from numpy import load,vstack,save,median

thisdir = os.getcwd()
os.chdir("/home/durant")
from pyraf import iraf
iraf.cd(thisdir)
iraf.digiphot()
iraf.daophot()
import pyraf
import pyfits
import numpy as n

def shift_file_coords(filename,xshift,yshift,output,sort=None):
    """Understands filetypes: 2-column ascii numbers, .mag, .als, .pst.
NB: shift means where each image is, relative to the original (not where
it should be moved to).
"""
    if not(sort):
        sort = 'num'
        if filename.find('.mag')>0: sort = 'mag'
        if filename.find('.als')>0: sort = 'als'
        if filename.find('.pst')>0: sort = 'pst'
    if not(sort=='num' or sort=='mag' or sort=='als' or sort=='pst'):
        raise ValueError('Unknown input filetype: %s'%filename)
    if sort=='num': # shift 2-column numeric ASCII table
        x,y = load(filename,usecols=[0,1],unpack=True)
        x += xshift
        y += yshift
        X = vstack((x,y))
        save(output,X.transpose())
        return
    if sort=='mag': #shift a .mag photometry file
        fred = open(filename)
        freda= open(output,'w')
        for line in fred:
            if line.split()[-1]=='\\' and len(line.split())==9 and line[0]!='#':
                x = float(line.split()[0]) + xshift
                y = float(line.split()[1]) + yshift
                line = "%-14.3f %-11.3f"%(x,y)+line[21:]
            freda.write(line)                
    if sort=='als': #shift a .als DAOphot photometry file
        fred = open(filename)
        freda= open(output,'w')
        for line in fred:
            if line.split()[-1]=='\\' and len(line.split())==8 and line[0]!='#':
                x = float(line.split()[1]) + xshift
                y = float(line.split()[2]) + yshift
                line = line[:9] + "%-10.3f %-10.3f"%(x,y) + line[29:]
            freda.write(line)
    if sort=='pst': #shift a PSF star list for DAOphot
        fred = open(filename)
        freda= open(output,'w')
        for line in fred:
            if line[0]!="#":
                x = float(line.split()[1]) + xshift
                y = float(line.split()[2]) + yshift
                line = line[:9] + "%-10.3f %-10.3f"%(x,y) + line[29:]
            freda.write(line)    
    fred.close()
    freda.close()

def recentre(image,refcoordfile):
    """Returns improved shift by centroiding
on the reference star using phot. This can be VERY
sensitive to the parameters in centerpars."""
    xin,yin = load(refcoordfile,unpack=True)
    try:
        iraf.phot(image,refcoordfile,'temp.mag',inter="no",calgorithm='centroid',
              mode='h',verify='no',update='no',verbose='no')
        xout,yout=iraf.pdump('temp.mag','xcen,ycen','yes',Stdout=1)[0].split()
    except:
        print "Recentring failed on", image
        return 0.,0.
    xout,yout = float(xout),float(yout)
    return xout-xin,yout-yin

vary_par = 1.
vary_max = 10
vary_min = 6
vary_fwhm= 0
    
def setaperture(image,refstar):
    """Measure the FWHM of the reference star unsing simple DAOphot editor
    and then set the photometry aperture to this number"""
    x,y = load(refstar,unpack=True)
    fred = open('tempaperfile','w')
    fred.write("%f %f 100 a\nq"%(x,y))
    fred.close()
    try:
        output=iraf.daoedit(image,icomm='tempaperfile',Stdout=1,Stderr=1)
    except:
        print "Aperture setting failed on",image
        return
    FWHM = float(output[3].split()[4])
    iraf.photpars.apertures = min(max(FWHM*vary_par,vary_min),vary_max)
    iraf.daopars.fitrad = min(max(FWHM*vary_par,vary_min),vary_max)
    global vary_fwhm
    vary_fwhm = FWHM
    print "FWHM: ", FWHM, " aperture: ",iraf.photpars.apertures

def apphot(image,coords,refstar=None,centre=False,vary=False):
    """Apperture photometry with centering based on a reference star.
NB: centre refers to shifting the coordinates by centroiding on the
reference star; recentering on the final phot depends on
centerpars.calgorithm ."""
    iraf.dele('temp.mag*')
    if centre:
        xsh,ysh = recentre(image,refstar)
        print "Fine centring: ", xsh,ysh
    else: #no recentreing by reference star (but could still have calgorithm!=none)
        xsh,ysh = 0,0
    if vary:
        setaperture(image,refstar)
    shift_file_coords(coords,xsh,ysh,'tempcoords')
    iraf.phot(image,'tempcoords','temp.mag2',inter="no",
                  mode='h',verify='no',update='no',verbose='no')
    out = iraf.pdump('temp.mag2','id,flux,msky,stdev','yes',Stdout=1)
    return out

def psfphot(image,coords,pststars,refstar,centre=True,vary=False):
    """PSF photometry. Centering is through phot on refstar.
Assume coords is a .als file for now. Recentering is always done
for the reference star, never for the targets."""
    iraf.dele('temp.mag*')
    iraf.dele('temp.psf.fits')
    iraf.dele('temp.als')
    if centre:
        xsh,ysh = recentre(image,refstar) 
        print "Fine Centring: ", xsh,ysh
    else: xsh,ysh = 0,0
    if vary:
        setaperture(image,refstar)
    shift_file_coords(coords,xsh,ysh,'tempcoords2',sort='als')
    shift_file_coords(pststars,xsh,ysh,'temppst2',sort='pst')
    iraf.phot(image,'tempcoords2','temp.mag2',inter="no",calgorithm='none',
                  mode='h',verify='no',update='no',verbose='no')
    iraf.psf(image,'temp.mag2','temppst2','temp.psf','temp.mag.pst','temp.mag.psg',
             inter='no',mode='h',verify='no',update='no',verbose='no')
    iraf.allstar(image,'temp.mag2','temp.psf','temp.als','temp.mag.arj',"default",
                 mode='h',verify='no',update='no',verbose='no')
    out = iraf.pdump('temp.als','id,mag,merr,msky','yes',Stdout=1)
    return out   

def simplepsfphot(image,coords,psf,refstar,centre=True,vary=False):
    """PSF photometry, with a given PSF file in psf used for every image"""
    iraf.dele('temp.mag*')
    iraf.dele('temp.als')
    iraf.dele('temp.sub.fits')
    if centre:
        xsh,ysh = recentre(image,refstar) 
        print "Fine Centring: ", xsh,ysh
    else: xsh,ysh = 0,0
    if vary:
        setaperture(image,refstar)
    shift_file_coords(coords,xsh,ysh,'tempcoords2',sort='als')
    iraf.phot(image,'tempcoords2','temp.mag2',inter="no",calgorithm='none',
                  mode='h',verify='no',update='no',verbose='no')
    iraf.allstar(image,'temp.mag2',psf,'temp.als','temp.mag.arj','temp.sub.fits',
                 mode='h',verify='no',update='no',verbose='no')
    out = iraf.pdump('temp.als','id,mag,merr,msky','yes',Stdout=1)
    return out   

def custom1(filename): # for NACO timing mode cubes - removes horizontal banding
    #iraf.imarith(filename,'-','dark','temp')
    iraf.imarith(filename,'/','flatK','temp')
    im = pyfits.getdata('temp.fits')
    med = median(im.transpose())
    out = ((im).transpose()-med).transpose()
    (pyfits.ImageHDU(out)).writeto("temp2.fits",clobber=True)
    iraf.imdel('temp')
    iraf.imcopy('temp2[1]','temp')

def get_id(starid,output='output'):
    """from the output of the photometry, grab the magnitudes and magerrs of starid"""
    mag = load(output,usecols=[4+starid*4])
    merr= load(output,usecols=[5+starid*4])
    return mag,merr

def run(filelist,coords,refstar,shifts=None,centre=False,psf=False,pststars=None,
        ids=None,dark=0,flat=1,timestamp="TIME",output='output',custom_process=None,
        vary=False):
    """If psf==True, must include all extra par files.
If PSF is a filename (.psf.fits), this profileis used to fit every image.
Timestamp can be either a file of times (same length as filelist), a header
keyword, or an array of times.
The input list can include [] notation for multiple extensions or sections
of each file (incompatible with header-based time-stamps).
custom_process(file) is a function taking a filename (possible including [x]
syntax) and places a processed image in temp.fits."""
    output = open(output,'w')
    x = load(coords,usecols=[1])
    numstars = len(x)
    myfiles = open(filelist).readlines()
    myfiles = [myfiles[i][:-1] for i in range(len(myfiles))]
    if timestamp.__class__ == numpy.ndarray: #--sort out times--
        times = 1 #times=1 means we know the times beforehand
    elif len(glob(timestamp))>0:
        timestamp = load(timestamp,usecols=[0])
        times=1
    else:
        times=0 #times=0 mean find the time from each image
    if type(shifts)==type(" "): #--sort out shifts--
        xshifts,yshifts = load(shifts,unpack=True)#filename give, assuming 2 columns
        xshifts,yshifts = -xshifts,-yshifts #these are in the opposite sense to coords from stack
    elif n.iterable(shifts): 
        xshifts=n.array(shifts[0])  #for shifts given as arrays/lists
        yshifts=n.array(shifts[1])
    else:
        print "No shifts" #assume all shifts are zero
        xshifts = n.zeros(len(myfiles))
        yshifts = n.zeros(len(myfiles))
    for i,thisfile in enumerate(myfiles): #run!
        print i,thisfile
        if times:  
            time = timestamp[i] #known time
        else:
            time = pyfits.getval(thisfile,timestamp) #FITS keyword
        try:
            iraf.dele('temp.fits')
            if custom_process: #arbitrary subroutine to process a file -> temp.fits
                custom_process(thisfile)
            else: #typical dark/bias subtract and flatfield
                iraf.imarith(thisfile,'-',dark,'temp')
                iraf.imarith('temp','/',flat,'temp')
            shift_file_coords(coords,xshifts[i],yshifts[i],'tempcoords') #apply coarse shifts
            shift_file_coords(refstar,xshifts[i],yshifts[i],'tempref',sort='num')
            if psf:
                if psf is True: #full PSF fit
                    shift_file_coords(pststars,xshifts[i],yshifts[i],'temppst')
                    out=psfphot('temp.fits','tempcoords','temppst','tempref',centre,vary)
                else: #DAOphot with known PSF
                    out=simplepsfphot('temp.fits','tempcoords',psf,'tempref',centre,vary)
            else: #aperture photometry
                out=apphot('temp.fits','tempcoords','tempref',centre,vary=vary)
            output.write("%s %s %s "%(thisfile,time,vary_fwhm))
            myids = n.array([int(out[i].split()[0]) for i in range(len(out))])
            for i in ids or range(numstars):
                try: #search for each requested ID
                    foundid = find(myids==i)[0]
                    output.write(out[foundid]+" ")
                except: #ID not found
                    output.write(" 0 0 0 0 ")
            output.write("\n")
        except KeyboardInterrupt: #exit on Ctrl-C
            break
        except pyraf.irafglobals.IrafError, err:
            print "IRAF error ",err,thisfile
            break
        except ValueError, err:
            print "Value error ",err,thisfile
            raise
    output.close()
    #iraf.dele('temp*')    
