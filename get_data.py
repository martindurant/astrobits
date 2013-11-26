tables = { '2MASS':"&-source=II/246/out&-out=RAJ2000&-out=DEJ2000&-out=Jmag&-out=e_Jmag&-out=Hmag&\
-out=e_Hmag&-out=Kmag&-out=e_Kmag&-out=Qflg&-out=Rflg&-out=Bflg&-out=Cflg&-out=Xflg&-out=Aflg",
    'USNOB':"&-source=I/284/out&-out=RAJ2000&-out=DEJ2000&-out=R1mag&-out=B1mag&-out=Imag",
    'GSC2':"&-source=I/305/out&out=RAJ2000&-out=DEJ2000&-out=Vmag&-out=Bmag",
    'UCAC':"&-source=I/289/out&-out=RA(ICRS)&-out=DE(ICRS)&-out=UCmag0",
    'SLOAN':"&-source=II/282/out&out=RAJ2000&-out=DEJ2000&-out=umag&-out=gmag&-out=rmag&-out=imag&-out=zmag",
        }

import urllib
urlbase = """http://vizier.u-strasbg.fr/viz-bin/asu-fits?-out.max=9999&-out.form=ascii%20999%27filled&-order=I&-c.eq=J2000&-c.u=arcmin&-c.geom=r&-out.add=_r&-sort=_r"""

def get(r,table,outfile,identifier=None,cenRA=None,cendec=None,verbose=False):
    """Call Vizier to retrive data for calibration, centred on cenRA/dec or
    "identifier" in Simbad, within r arcminutes.
    "table" is the data to retrieve (2MASS|UCAC|GSC2|USNOB|SLOAN)"""
    try:
        querie = tables[table]
    except KeyError:
        raise KeyError("Tables implementes: 2MASS|UCAC|GSC2|USNOB|SLOAN")
    if not(identifier is None):
        coord = '&-c=%s'%identifier.replace(" ","+") + "&-c.r=%f"%r
    else:
        coord = "&-c=%f%+f"%(cenRA,cendec) + "&-c.r=%f"%r
    url = urlbase+querie+coord
    if verbose: print url
    urllib.urlretrieve(url,outfile)

import xml.dom.minidom
from xml.dom.minidom import Node

def resolve(name):
    """Get coordinates for 'name' from simbad."""
    name2 = name.replace('+',u'%2B')
    name3 = name2.replace(' ','+')
    url = r"http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oxp/SNV?%s"%name3
    try:
        out,msg = urllib.urlretrieve(url)
        doc = xml.dom.minidom.parse(out)
        r = float(doc.getElementsByTagName('Resolver')[0].getElementsByTagName('jradeg')[0].childNodes[0].data)
        d = float(doc.getElementsByTagName('Resolver')[0].getElementsByTagName('jdedeg')[0].childNodes[0].data)
        return r,d
    except:
        print "Fetch failed",url

