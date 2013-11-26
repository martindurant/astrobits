import numpy as n

from pyraf import iraf
iraf.sdsdas()
iraf.hst_calib()
iraf.synphot()

calcphot V "rn(icat(k93models,i,0.0,3.0)*ebmv($0),band(wfc3,uvis2,f225w,cal),184,counts)"  fnu  out=b0656_aqu_test.fits wave=wave_1300_20000_5.fits vzero="0.02-0.2x0.02"  >> test.dat

def get_magni
