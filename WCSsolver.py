"""Given a FITS table of sources around a position, find an
astrometric solution for an image in an interactive way. The
image must have an associated .als photometry file"""

import CMD
import pyfits
from pyraf import iraf
iraf.martin()
from pylab import clf,plot,figure,imshow

fields = {"GSC":("RA(ICRS)","DE(ICRS)")}

def findsolution(image,alsfile,datafile,source="GSC"):
    stars = CMD.star_photometry(alsfile,image,False)
    datfile = pyfits.open(datafile)
    mydata=datfile[-1].data
    datfile.close()
    iraf.display(image,1)
    temp=iraf.pdump(alsfile,'xcen,ycen,id','yes',Stdout=1)
    iraf.tvmark(1,"STDIN",Stdin=temp)
    rafield,decfield = fields[source]
    ra = mydata.field(rafield)
    dec= mydata.field(decfield)
    raint = int(min(ra))
    decint= int(min(dec))
    print "ra offset by %i, dec by %i\n"%(raint,decint)
    clf()
    plot(ra-raint,dec-decint,'k+')
    print "Require 4 stars matched to make simple solution"
    print "Enter co-ordinates as acurately as you can:-"
    fred = open('ccmapfile','w')
    for i in range(4):
        idin=input("Star %i ID:  "%i)
        mystar=stars[idin]
        rain=input("Star %i RA:  " %i) + raint
        decin=input("Star %i DEC: "%i) + decint
        radiff = ra-rain
        decdiff = dec-decin
        radecdiff = radiff**2+decdiff**2
        locmin = radecdiff.argmin()
        raout = ra[locmin]
        decout = dec[locmin]
        fred.write("%7.2f %7.2f %10.6f %10.6f\n"%(mystar[1],mystar[2],raout,decout))
    fred.close()
    iraf.ccmap('ccmapfile','database',images=image,xcolumn=1,ycolumn=2,lngcolumn=3,
               latcolumn=4,lngunits='degrees',latunits='degrees',insystem='j2000',
               fitgeometry='rxyscale')
    findbettersolution(image,alsfile,datafile,source=source)
    
def findbettersolution(image,alsfile,datafile,source="GSC",maxmag=30,minmag=10,chicut=3):
    stars = CMD.star_photometry(alsfile,image,True)
    for star in stars.data:
        if star[3]>maxmag or star[3]<minmag or star[6]>chicut:
            star[7]=-5000 #don't find stars outside mag range
    datfile = pyfits.open(datafile)
    mydata = datfile[-1].data
    datfile.close()
    rafield,decfield = fields[source]
    fred=open('ccmapfile','w')
    for refstar in mydata:
        ra = refstar.field(fields[source][0])
        dec= refstar.field(fields[source][1])
        mystar = stars.starbyWCS(ra,dec,limit=2./3600)
        if type(mystar) is not type(None):
            fred.write("%7.2f %7.2f %10.6f %10.6f\n"%(mystar[1],mystar[2],ra,dec))
    fred.close()
    iraf.ccmap('ccmapfile','database',images=image,xcolumn=1,ycolumn=2,lngcolumn=3,
               latcolumn=4,lngunits='degrees',latunits='degrees',insystem='j2000',
               fitgeometry='general')
    
