import os,string,sys
import pylab
import numpy as numarray
from numpy import zeros,where,add,sqrt,array,any,nan,arange
thisdir = os.getcwd()
os.chdir('/home/durant')
from pyraf import iraf
os.chdir(thisdir)
iraf.cd(thisdir)
iraf.stsdas()
iraf.digiphot()
iraf.daophot()
import glob
import pyfits
import astLib


class star_photometry(object):
    """reads photometry information from IRAF .als or .mag file, and
    returns information from it as arrays of magnitudes etc, including
    conversions of co-ordinates etc.
    self.data is a float array of the data, containing:
    id,x,y,mag,magerr,chi,ra,dec.
    self.image is the corresponding FITS file with WCS information
    self.photfile is the original file
    self.nstar is number of stars (=len(self.data))"""
    mynames={"id":0,"x":1,"y":2,"mag":3,"magerr":4,
             "sky":5,"chi":6,"ra":7,"dec":8} #record structure
    def __init__(self,photfile,imagefile=None,calcradec=True):
        """load data from photfile, and store in self.stars,
        associate with imagefile if it exists, or try to find it"""
        self.image=None
        try:
            photfile=glob.glob(photfile)[0]
        except:
            raise IOError("File not found: %s"%photfile)
        self.photfile=photfile
        try:            
            rootname=photfile.split('.')[0]
            if imagefile:
                imagefile=glob.glob(imagefile)[0]
            else:
                imagefile=(glob.glob(rootname+'.fits') or glob.glob(rootname+"*"+"fits"))[0]
                print "Setting image file to: %s"%imagefile
            self.image=imagefile
        except:
            print "No image file found - set self.image manually"
        fred=open(photfile)
        mylines=fred.readlines()
        fred.close()
        if 'als' in photfile.split("."):
            nstar=(len(mylines)-44)//2 #ignore comment lines, two lines per star
            self.nstar=nstar
            self.data=zeros(shape=(nstar,9))*1.0
            for index in range(nstar):
                line=mylines[index*2 + 44]
                words=line.split()
                self.data[index,:6]=[float(words[i]) for i in range(6)]
                line=mylines[index*2 + 45]
                words=line.split()
                self.data[index,6]=float(words[1])
        elif 'mag' in photfile.split("."):
            nstar=(len(mylines)-75)//5 #ignore comment lines, two lines per star
            self.nstar=nstar
            self.data=zeros(shape=(nstar,9))*1.0
            for index in range(nstar):
                line=mylines[index*5 + 75]
                self.data[index,0] = int(line.split()[3])
                words=mylines[index*5 + 76].split()
                self.data[index,1:3] = (float(words[0]),float(words[1]))
                line=mylines[index*5 + 77]
                self.data[index,5] = float(line.split()[0])
                words = mylines[index*5 + 79].split()
                try:
                    self.data[index,3:5]= (float(words[4]),float(words[5]))
                except ValueError:
                    self.data[index,3:5]=(nan,nan)
        else:
            raise IOError("File %s of uknown photometry type\ntry .mag or .als"%photfile)
        if calcradec and self.image:
            self.calcradec()
        
    def __repr__(self):
        if self.image:
            return "Photometry from file %s and image %s: %i stars"%(
                self.photfile,self.image,self.nstar)
        else:
            return "Photometry from file %s, no image set: %i stars"%(
                self.photfile,self.nstar)

    def allcal(self,box,table,identifier=None,cenRA=None,cendec=None):
        """Call Vizier to retrive data for calibration, centred on cenRA/dec or
        "identifier" in Simbad, of box arcminutes wide.
        "table" is the data to retrieve (2mass|ucac|GSC2|USNOB)"""
        import get_data
        if cenRA is None or cendec is None:
            cenRA = pyfits.getval(self.image,"RA")
            cendec= pyfits.getval(self.image,"DEC")
        rafield,decfield,mags=getdata.getfile(box,table,'allcal.fits',identifier,cenRA,cendec)
        self.autocal('allcal.fits',rafield,decfield)
        self.ast_by_file('allcal.fits',rafield,decfield)
        mymag = raw_input('Chose mag from %s"%mags: ')
        if len(mymag)>0:
            out=selt.cal_by_file('allcal.fits',mag,rafield,decfield)
            pylab.hist(out,bins = arange(median(out)-1,median(out)+1,0.05))
            return out

    def photautocal(self,reference,reflimit=30,selflimit=30):
        """Find the rough WCS for this image, based on reference star_phootmetry,
        by cross-identifying stars by ID in DS9"""
        iraf.display(self.image,1)
        temp = iraf.pdump(self.photfile,"xcen,ycen,id","MAG<%f"%selflimit,Stdout=1)
        iraf.tvmark(1,"STDIN",Stdin=temp)
        iraf.display(reference.image,2)
        temp = iraf.pdump(reference.photfile,"xcen,ycen,id","MAG<%f"%reflimit,Stdout=1)
        iraf.tvmark(2,"STDIN",Stdin=temp)
        fred = open('tempautocal','w')
        input1=None
        input2=None
        print "Cross identify at least five stars, <q> to continue"
        while (input1!='q')and(input2!='q'):
            print
            input1 = raw_input("Star ID on Frame 1: ")
            input2 = raw_input("Matches star ID on Frame 2: ")
            try:
                star1 = self[int(input1)]
                star2 = reference[int(input2)]
            except:
                print "Known star integers required"
                continue
            fred.write("%f %f %f %f\n"%(star1[1],star1[2],star2[7],star2[8]))
        fred.close()
        iraf.ccmap(input='tempautocal',images=self.image,lngunits='degrees',
                   database='database',fit='rxyscale',update='yes',inter='yes')
        print "Recalculating data"
        self.calcradec()
        print "If happy, now run ast_by_photometry"

    def ast_by_photometry(self,phot_data,mag_low=15,mag_high=20,limit=1./3600):
        """Take another star_photometry and refine the WCS in this image
        by matching stars. photautocal can be used to find the rough
        WCS"""
        if phot_data.image==None:
            print "Need calibrated photometry to calibrate from"
        RA = phot_data['ra']
        dec= phot_data['dec']
        mag= phot_data['mag']
        ind= (mag>mag_low)*(mag<mag_high)
        self.ast(RA[ind],dec[ind],limit)

    def ast(self,RA,dec,limit):
        """Does the astrometric calibration. You can call this directly, if you
        want to load data and filter before calibrating"""
        fred = open('tempcoordfile','w')
        for i in range(len(RA)):
            star = self.starbyWCS(RA[i],dec[i],limit=limit)
            if not(star is None):
                fred.write("%f %f %f %f\n"%(star[1],star[2],RA[i],dec[i]))
        fred.close()
        iraf.ccmap(input='tempcoordfile',images=self.image,lngunits='degrees',
                   database='database',fit='rxyscale',update='yes',inter='yes')
        print "Recalculating data"
        self.calcradec()

    def cal_by_photometry(self,phot_data,mag_low,mag_high,limit=1./3600):
        if phot_data.image==None:
            print "Need calibrated photometry to calibrate from"
        RA = phot_data['ra']
        dec= phot_data['dec']
        mag= phot_data['mag']
        ind= (mag>mag_low)*(mag<mag_high)
        return self.cal(RA[ind],dec[ind],mag[ind],limit)

    def cal(self,RA,dec,mag,limit):
        """Does the photometric calibration. You can call this directly, if you
        want to load data and filter before calibrating"""
        results = []
        results2= []
        for i in range(len(RA)):
            star = self.starbyWCS(RA[i],dec[i],limit=limit)
            if not(star is None):  #this loop is very slow for large numbers!
                results.append(star[3]-mag[i])
                results2.append(star[3])
        return array(results),array(results2)
        
    def autoast(self,infile='asu.fit',RAfield='RAJ2000',decfield='DEJ2000',magfield='Vmag',magoff=None,
                markmax=30,crosses=False):
        """Calibrate image from file of data (2mass etc) by making initial astrometry
        by selecting a handful of stars, them using the other cal routines
        RAfield: column with RA (decimal!).
        decfield: column with dec (decimal!).
        magfield: column to calibrate magnitude by"""
        global xpoint,ypoint
        if self.image is None:
            print "No image to calibrate data by"
            return
        RA,dec,mag = load_data(infile,RAfield,decfield,magfield)
        fig = pylab.figure(figsize=(8,8))
        pylab.clf()
        try:
            magoff = magoff or mag.max()
            pylab.scatter(RA[mag>0],dec[mag>0],s=10**(0.4*(magoff-mag[mag>0])),hold=0,alpha=0.5)
        except:
            crosses=True
        if crosses: pylab.plot(RA,dec,'k+')
        pylab.axis([RA.max(),RA.min(),dec.min(),dec.max()])
        pylab.draw()
        cid = fig.canvas.mpl_connect('button_press_event',self._on_press)
        iraf.display(self.image,1)
        iraf.tvmark(1,"STDIN",Stdin=iraf.pdump(self.photfile,"xcen,ycen,id","MAG<%f"%markmax,Stdout=1))
        print "Pick five stars with middle button, then enter their IDs in order,\n<q> to continue"
        ids = []; self.xpoint=[]; self.ypoint=[]
        select = None
        while not(select=='q'):        
            select=raw_input('Star ID: ')
            try: ids.append(int(select))
            except: pass
        fig.canvas.mpl_disconnect('button_press_event')
        fred = open('autocalcoords','w')
        for i in range(min(len(ids),len(self.xpoint))):
            try:
                x,y = self[ids[i]][[1,2]]
                locmin = ((self.xpoint[i]-RA)**2+(self.ypoint[i]-dec)**2).argmin()
                fred.write('%f %f %15.12e %15.12e\n'%(x,y,RA[locmin],dec[locmin]))
            except:
                print "Coord pair %i failed"%i
        fred.close()
        iraf.ccmap(input='autocalcoords',images=self.image,lngunits='degrees',
                   database='database',fit='rxyscale',update='yes',inter='yes')
        print "Recalculating..."
        self.calcradec()
        print "If happy, now run ast_by_file"        
    
    def _on_press(self,event):
        """Helper function for interactive clicking in autocal"""
        if event.button==2 and event.xdata!=0:
            self.xpoint.append(event.xdata)
            self.ypoint.append(event.ydata)

    def cal_by_file(self,tablefile='asu.fit',magfield='Vmag',rafield='RAJ2000',decfield='DEJ2000',
                    mag_low=10,mag_high=20,limit=.1/3600,clip=True):
        """For an input FITS with columns RA, dec, mag; find closest star for our data by
        WCS for each star and calculate (but do not set) magnitude offset. magfield is
        the column to calibrate on (e.g., Kmag in 2mass). Rejects reference magnitudes
        outside of maglow/high.
        Returns list of mag offsets, so you can make your own cuts/average."""
        RA,dec,mag = load_data(tablefile,rafield,decfield,magfield)
        ind= (mag>mag_low)*(mag<mag_high)
        if clip:
            ind *= (RA>self['ra'].min())*(RA<self['ra'].max())*(
                dec>self['dec'].min())*(dec<self['dec'].max())
        return self.cal(RA[ind],dec[ind],mag[ind],limit)

    def ast_by_file(self,tablefile='asu.fit',rafield='RAJ2000',decfield='DEJ2000',limit=.1/3600,clip=True):
        """Given file of RA, dec (2mass etc); cross-identify stars,
        call ccmap and then improve the fit using all the stars in the table
        centroiding. 
        Image needs rough WCS to proceed, e.g. from autoast."""
        data = pyfits.getdata(tablefile)
        RA = data.field(rafield)
        dec = data.field(decfield)
        if clip:
            ind = (RA>self['ra'].min())*(RA<self['ra'].max())*(
                dec>self['dec'].min())*(dec<self['dec'].max())
        self.ast(RA,dec,limit)

    def copy(self):
        """make a new version (for masking etc)"""
        newan = star_photometry(self.photfile,self.image,calcradec=0)
        newan.data = self.data.copy()
        return newan

    def xyoffset(self,x,y):
        """shift all co-ordinates by x,y"""
        self.data[:,1]+=x
        self.data[:,2]+=y

    def magoffset(self,mag):
        """add magnitude zero point"""
        self.data[:,3]+=mag

    def calcradec(self):
        """calculate ra and dec in hours/degrees from x,y coords
        in data, image WCS and xy2ra from StSci"""
        if self.image is None:
            print "No photometry image set, not doing WCS"
            return
        self.wcs = astLib.astWCS.WCS(self.image)
        for star in self.data:
            star[7],star[8] = self.wcs.pix2wcs(star[1],star[2])

    def starbyxy(self,x,y,limit=10):
        """for coord(s) x,y, return the nearest star(s)"""
        xdiff = self.data[:,1]-x
        ydiff = self.data[:,2]-y
        xydiff = xdiff**2+ydiff**2
        locmin = xydiff.argmin()
        if sqrt(xydiff[locmin])>limit:
            return None
        return self.data[locmin]

    def starbyWCS(self,ra,dec,limit=.1/3600.):
        """for coord(s) ra,dec in decimal degrees, return
        the nearest star(s). Distance limit set equal to 0.1"."""
        radiff = self.data[:,7]-ra
        decdiff = self.data[:,8]-dec
        radecdiff = radiff**2+decdiff**2
        locmin = radecdiff.argmin()
        if sqrt(radecdiff[locmin])>limit:
            return None
        return self.data[locmin]

    def __getitem__(self,starid):
        """return star record by id"""
        if type(starid) is type(" "):
            return self.data[:,self.mynames[starid]]
        if type(starid) is type(1):
            return self.data[where(self.data[:,0]==starid)][0]
        start=starid.start or 0
        if start<0: start += len(self.data)
        stop =starid.stop or max(self.data[0])
        if stop<0: stop += len(self.data)
        if starid.step:
            ids = range(starid.start,starid.stop,starid.step)
        else:
            ids = range(starid.start,starid.stop)
        out=[]
        for myid in ids:
            if len(self[myid])>0:
                out.append(self[myid])
        return out

def load_data(tablefile='asu.fits',rafield='RAJ2000',decfield='DEJ2000',magfield='Vmag'):
    data = pyfits.getdata(tablefile)
    RA = data.field(rafield)
    dec = data.field(decfield)
    try:
        mag = data.field(magfield)
    except:
        mag = None
    return RA,dec,mag
    
 
def CMD(refstars,findstars,limit=3,chicut=2,method='xy',verbose=False,idout=False):
    """take a set of stars, and return the magnitudes of the nearest stars
    in the other set.
    refstars: the stars to search for (usually the smaller number)
    findstars: the stars to search from
    limit: how far (in x,y) appart stars are allowed to be
    chicut: ignore all stars with chi bigger than this
    idout: whether to include the column of IDs in refstars
    returns mag1,merr1,mag2,merr2 for found stars"""
    mag1=[]; merr1=[]
    mag2=[]; merr2=[]
    ids=[]
    for mystar in refstars.data:
        if method=='xy':
            foundstar = findstars.starbyxy(mystar[1],mystar[2],limit=limit)
        else:
            foundstar = findstars.starbyWCS(mystar[7],mystar[8],limit=limit)
        if not(foundstar is None):
            if mystar[6]>chicut or foundstar[6]>chicut: continue
            mag1.append(mystar[3])
            merr1.append(mystar[4])
            mag2.append(foundstar[3])
            merr2.append(foundstar[4])
            ids.append(mystar[0])
            if verbose:
                print "Star %i matches Star %i"%(mystar[0],foundstar[0])
    if idout:
        return array(mag1),array(merr1),array(mag2),array(merr2),array(ids)
    return array(mag1),array(merr1),array(mag2),array(merr2)

def CCD(refstars,findstars1,findstars2,limit=3,chicut=2,method='xy',verbose=False,idout=False):
    """take a set of stars and return the magnitudes of the nearest stars
    in both other sets.
    refstars: the stars to search for (usually the smaller number)
    findstars1/2: the stars to search from
    limit: how far (in x,y) appart stars are allowed to be
    method: x/y for pixel co-ords, WCS otherwise
    chicut: ignore all stars with chi bigger than this
    idout: whether to include the column of IDs in refstars
    returns mag1,merr1,mag2,merr2,mag3,merr3 for found stars"""
    mag1=[]; merr1=[]
    mag2=[]; merr2=[]
    mag3=[]; merr3=[]
    ids=[]
    for mystar in refstars.data:
        if method=='xy':
            foundstar1 = findstars1.starbyxy(mystar[1],mystar[2],limit=limit)
            foundstar2 = findstars2.starbyxy(mystar[1],mystar[2],limit=limit)
        else:
            foundstar1 = findstars1.starbyWCS(mystar[7],mystar[8],limit=limit)
            foundstar2 = findstars2.starbyWCS(mystar[7],mystar[8],limit=limit)
        if not(foundstar1 is None) and not(foundstar2 is None):
            if mystar[6]>chicut or foundstar1[6]>chicut or foundstar2[6]>chicut: continue
            mag1.append(mystar[3])
            merr1.append(mystar[4])
            mag2.append(foundstar1[3])
            merr2.append(foundstar1[4])
            mag3.append(foundstar2[3])
            merr3.append(foundstar2[4])
            ids.append(mystar[0])
            if verbose:
                print "Star %i matches Star %i and Star %i"%(
                mystar[0],foundstar1[0],foundstar2[0])
    if idout:
        return array(mag1),array(merr1),array(mag2),array(merr2),array(mag3),array(merr3),array(ids)
    return array(mag1),array(merr1),array(mag2),array(merr2),array(mag3),array(merr3)

def catalogue(inputlist,limit=.1/3600):
    """Take a list of photometry objects and make a catalogue by cross-matching the
    known RA and dec within "limit" distance. Outputs RA,DEC,(m1,e1...) In the first
    version, this will produce one line per star in the first list"""
    numphot = len(inputlist)
    if numphot<2 or type(inputlist)!=type((1,)):
        print "Input >2 photometry objects as a list or tuple"
    numstars = len(inputlist[0].data)
    ra = inputlist[0].data[:,7]; dec = inputlist[0].data[:,8]
    output = zeros((numphot,2,numstars))*1.0
    for i in range(numstars):
        output[0][0][i] = inputlist[0].data[i][3]
        output[0][1][i] = inputlist[0].data[i][4]
        for j in arange(numphot-1)+1:
            star = inputlist[j].starbyWCS(ra[i],dec[i],limit=limit)
            if not(star is None):
                output[j][0][i] = star[3]
                output[j][1][i] = star[4]
    return ra,dec,output

