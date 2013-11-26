import Numeric as N
import numarray as num
from Scientific.Functions.LeastSquares import leastSquaresFit as LS
from cashfitter import CashFit as CF
from pylab import *
import xspecs
from Scientific.Statistics import standardDeviation
import Kedge_cash

__all__ = ['edge','deextinct','plot_fe','plot_ne','plot_mg','plot_si','get_edge',
           'auto_edge','get_edge','MC1','MC2']

# LeastSquaresFit assumes funtions defined in Numeric

cash = 0

data = {
    'C':(10.5274  ,297.37),
    'C+':(9.0638  ,315.15),
    'N':(7.5208 ,412.36),
    'N+':(6.6628  ,432.2),
    'O':(5.6419    ,546.02),
    'O+':(5.1187  ,567.94),
    'Ne':(3.5231  ,869.4),
    'Ne+':(3.265  ,895.45),
    'Mg':(2.191   ,1310.58),
    'Mg+':(2.1411,1322.29),
    'Si':(1.4759   ,1848.62),
    'Si+':(1.4481 ,1862.2),
    'FeL':(4.936 ,730.), #for iron L-shell edge (p electrons only)
    'FeK':(0.3331,7123.6)
         }  #sigma_K(Ethres) is in units of 10^-19 cm2 (first datum in each row)

spectrum=None
l_edge=None
sig=None
edge_power=2.5 # cross-section falls off as (l/l_e)**edge_power past edge

def get_edge(l,flux,err,edge,lowl,highl):
    """Fits a photoelectric edge to the section of the l/flux spectrum
    falling between lowl and highl.
    l: wavelength in angstrom
    flux: some measure of flux (erg/s/cm2/A?)
    edge: one of the edges in data (O,Ne,FeL,Si...)"""
    global l_edge,sig,sigma_K_Et,E_K
    sigma_K_Et,E_K = data[edge]
    l_edge = 2.998e18/(E_K *1.602e-19/6.626e-34)
    l=num.array(l)
    indeces=(l>lowl)*(l<highl)
    if edge=='FeL':
        l_edge=17.52
        indeces *= (l>17.2)*(l<17.56)-1
    elif edge=='Ne':
        l_edge=14.31    #http://space.mit.edu/HETG/Reports/HETG_Report_SciFeb03.html
        indeces *= (l<14.44)*(l>14.66)-1
        indeces *= (l>13.4)*(l<13.5)-1
    elif edge=='Mg':
        l_edge=9.5
    elif edge=='Si':
        l_edge=6.72
        indeces *= (l<6.73)*(l>6.61)-1
    elif edge=="O":
        l_edge=23.1
        indeces *= (l<23.6)*(l>23.25)-1
        indeces *= (l<23.1)*(l>22.5)-1
    lne = l[num.where(indeces)]
    fne = flux[num.where(indeces)]
    ene = err[num.where(indeces)]
    sig = sigma_K_Et*1e-19
    return lne,fne,ene

def fit_edge(l,flux,err,guess_norm=None, guess_alpha=None, guess_N=None):
    """Does edge fitting. Either called graphically from edge() one edge at
    a time or multiple times from auto_edge()"""
    inputs = N.zeros((len(l),3),typecode='d')
    inputs[:,2]=err
    inputs[:,0]=l
    inputs[:,1]=flux
    guess_N    =   guess_N or 0.1/sigma_K_Et
    guess_alpha = guess_alpha or -1.
    guess_norm = guess_norm or N.average(flux)/E_K**guess_alpha
    if cash==0:
        return LS(func_for_fit,(guess_alpha,guess_N,guess_norm),inputs)
    else:
        return Kedge_cash.cashfit(l,flux,err,sigma_K_Et,l_edge)
#output[0][0]=>alpha; output[[0][1]=>column N; output[0][2]=>norm A;  output[1]=>chi2

def auto_edge(l,flux,err,edge,lowl,highl):
    """Gets value of column for 'edge' element, without the graphical output.
    'results' is available globally (Kedge.results) if needed"""
    global results
    l,flux,err = get_edge(l,flux,err,edge,lowl,highl)
    results = fit_edge(l,flux,err)
    return results[0][1]

def edge(l,flux,err,edge,lowl,highl):
    """ Takes a small spectral range (E,flux,err), and fits a column density of
    element 'edge' to it, assuming the spectrum is a powerlaw in this region.
    'edge' is one of the values in the data dict above. This is the graphical
    version. For rapid calculation, use auto_edge()"""
    figure(1)
    clf()
    try:axis([min(l)-0.05,max(l)+0.05,min(flux)*0.9,max(flux)*1.1])
    except: pass
    l,flux,err = get_edge(l,flux,err,edge,lowl,highl)
    loglog(l[where(err<9)],flux[where(err<9)],'ko',markersize=4)
    #loglog(l[where(err>9)],flux[where(err>9)],'ko',mfc='w',markersize=4)
    print "Numer of data points:", len(l)
    results = fit_edge(l,flux,err)
    specl = arange(min(l),max(l),(max(l)-min(l))/10000) #for plotting only
    spectrum = results[0][2] * (specl/l_edge)**results[0][0] * N.where(specl>l_edge,1,\
            N.exp(-(specl/l_edge)**3*sig*results[0][1]))
    spectrum2 = results[0][2] * (specl/l_edge)**results[0][0]
    plot(specl,spectrum,'r-',specl,spectrum2,'r:')
    figure(2)
    clf()
    chis = chi_err(l,flux,err,results[0],arange(-1e17,3e18,1e16))
    plot(arange(-1e17,3e18,1e16),chis)
    return results

def func_for_fit(para,l):
    """Fitting function, with p = (power law index, column density, normalisation)"""
    alpha = para[0]
    column = para[1]
    A = para[2]
    if l>l_edge: factor=1
    else: factor=N.exp(-(l/l_edge)**edge_power*sig*column)
    return A * (l/l_edge)**alpha * factor

def MC1(l,f,err,edge,lowl,highl,reps):
    """Bootstrap method for producing a distribution of columns about the best
    fit. reps is the number of runs to make"""
    global output
    output = zeros(reps)*1.0
    index = zeros(reps)*1.0 #if you want to store both parameters
    npoints= len(l)
    for count in range(len(output)):
        if count % 100 == 0: print count
        indeces = (rand(npoints) * npoints).astype(int16)
        lsynth = l[indeces]
        fsynth = f[indeces]
        esynth = err[indeces]
        output[count] = auto_edge(lsynth,fsynth,esynth,edge,lowl,highl)
        index[count] = results[0][0] #you want indeces?
    return output,index #you want indeces?

def MC2(l,f,err,edge,lowl,highl,reps):
    """Monte-carlo method, creates datasets by applying the variance in the
    data to the best-fit model. reps is the number of runs to make"""
    global output
    auto_edge(l,f,err,edge,lowl,highl) #get best-fit
    lne,fne,ene = get_edge(l,f,err,edge,lowl,highl)
    spectrum = results[0][2] * (lne/l_edge)**results[0][0] * N.where(lne>l_edge,1,\
            N.exp(-(lne/l_edge)**3*sig*results[0][1]))
    goodpoints = where(ene<9)
    sd = standardDeviation(fne[goodpoints] - spectrum[goodpoints])
    print "Standard deviation of residuals:",sd
    output = zeros(reps)*1.0
    indeces = zeros(reps)*1.0 #if you want to store both parameters
    for count in range(reps):
        if count%100 == 0 : print count
        fnew = spectrum + randn(len(spectrum))*sd
        output[count] = auto_edge(lne,fnew,ene,edge,lowl,highl)
        indeces[count] = results[0][0] #you want indeces?
    return output,indeces #you want indeces?
                               

def chi_err(l,f,err,fit_pars,columns):
    """Produces a nice plot of chi^2 versus different values of the column for
    estimation of the 1-sigma error. l,f and err are the spectral data, fit pars
    are from fit_edge above and columns is the range of values to calculate chi^2
    for"""
    alpha,column,norm = fit_pars[0],fit_pars[1],fit_pars[2]
    alt_l = N.where(l>l_edge,1,l)
    factor = N.exp(N.multiply.outer(columns*sig,-(alt_l/l_edge)**edge_power))
    spectra = norm * (l/l_edge)**alpha * factor
    deviations = (spectra-f)**2/err**2
    chis = sum(N.transpose(deviations))
    return chis

def wabs(f,Nh,l=None,e=None):
    """
    Apply Morrison&McCammon, 1982 to wavelengths l(A) or e(keV)
    to de-extinct spectrum f. Use negative Nh to extinct model
    spectra"""
    mult = f/f
    if type(l)<>type(None):
        e = (1.986e-15/l) / 1.602e-16
    if type(e)==type(None):
        raise TypeError("...and the spectral range?")
    index=where((e>=0.03)*(e<0.1))
    mult[index] = (1e-24 * (17.3 + 608.1*e - 2150*e**2)/e**3)[index]
    index=where((e>=0.1)*(e<.284))
    mult[index] = (1e-24 * (34.6 + 267.9*e - 476.1*e**2)/e**3)[index]
    index=where((e>=0.284)*(e<.4))
    mult[index] = (1e-24 * (78.1 + 18.8*e + 4.3*e**2)/e**3)[index]
    index=where((e>=0.4)*(e<.532))
    mult[index] = (1e-24 * (71.4 + 66.8*e - 51.4*e**2)/e**3)[index]
    index=where((e>=0.532)*(e<.707))
    mult[index] = (1e-24 * (95.5 + 145.8*e - 61.1*e**2)/e**3)[index]
    index=where((e>=0.707)*(e<.867))
    mult[index] = (1e-24 * (308.9 - 380.6*e + 294*e**2)/e**3)[index]
    index=where((e>=0.867)*(e<1.303))
    mult[index] = (1e-24 * (120.6 + 169.3*e - 47.7*e**2)/e**3)[index]
    index=where((e>=1.303)*(e<1.84))
    mult[index] = (1e-24 * (141.3 + 146.8*e - 31.5*e**2)/e**3)[index]
    index=where((e>=1.84)*(e<2.471))
    mult[index] = (1e-24 * (202.7 + 104.7*e - 17*e**2)/e**3)[index]
    index=where((e>=2.471)*(e<3.21))
    mult[index] = (1e-24 * (342.7 + 18.7*e - 0*e**2)/e**3)[index]
    index=where((e>=3.21)*(e<4.038))
    mult[index] = (1e-24 * (352.2 + 18.7*e - 0*e**2)/e**3)[index]
    index=where((e>=4.038)*(e<7.111))
    mult[index] = (1e-24 * (433.9 - 2.4*e - 0.75*e**2)/e**3)[index]
    index=where((e>=7.111)*(e<8.331))
    mult[index] = (1e-24 * (629 + 30.9*e - 0*e**2)/e**3)[index]
    index=where((e>=8.331)*(e<10))
    mult[index] = (1e-24 * (701.2 + 25.2*e - 0*e**2)/e**3)[index]
    return f * exp(mult*Nh)

def deextinct(l,f,ocolumn=None,fecolumn=None,necolumn=None
              ,mgcolumn=None,sicolumn=None,Nh=0,normdepth=None):
    """take the values for column densities (cm-2) and de-extinct
    the lambda,flux spectrum using them. normdepth stands for all
    the optical depth (*lambda**3) below 25A, normalised at 25A. Uses
    abundances from Asplund et al (2004), for unspecified columns.
    Use negative columns to extinct model spectra."""
    global depth
    normdepth = normdepth or Nh*7.6e-22
    ocolumn = ocolumn or Nh*4.57e-4
    fecolumn = fecolumn or Nh*2.81e-5
    necolumn = necolumn or Nh*6.92e-5
    mgcolumn = mgcolumn or Nh*3.39e-5
    sicolumn = sicolumn or Nh*3.24e-5
##    ocolumn = ocolumn or Nh*7.41e-4  # old abundances as in Morrison&McCammon
##    fecolumn = fecolumn or Nh*3.31e-5
##    necolumn = necolumn or Nh*1.38e-4
##    mgcolumn = mgcolumn or Nh*3.98e-5
##    sicolumn = sicolumn or Nh*3.72e-5
    depth = l/l
    depth = exp((l/25)**3*normdepth)
    depth[where(l<=23.1)] *= exp((l[where(l<=23.1)]/23.1)**2.5*ocolumn*5.6419e-19)
    depth[where(l<=17.52)] *= exp((l[where(l<=17.52)]/17.52)**2.5*fecolumn*4.936e-19) #L-edges
    depth[where(l<=14.31)] *= exp((l[where(l<=14.31)]/14.31)**2.5*necolumn*3.5231e-19)
    depth[where(l<=9.5)]   *= exp((l[where(l<=9.5)]/9.5)**2.5*mgcolumn*2.191e-19)
    depth[where(l<=6.72)]  *= exp((l[where(l<=6.72)]/6.72)**2.5*sicolumn*1.4759e-19)
    depth[where(l<=1.74)]  *= exp((l[where(l<=1.74)]/1.74)**2.5*fecolumn*0.3331e-19) #K-edge
    f_fixed = f*depth
    return f_fixed

def plot_fe():
    os.chdir('/data/xmm/0142')
    l,f,e=xspecs.data_in('0142_1.fits')
    close()
    figure(1,(7,7))
    subplot(211)
    edge(l,f,e,"FeL",16,19)
    ylim(0.7e-4,0.7e-3)
    xticks(arange(16,20))
    os.chdir('../2259')
    l,f,e=xspecs.data_in('2259_1.fits')
    figtext(0.6,0.8,"4U~0142+61",size=18)
    subplot(212)
    edge(l,f,e,"FeL",16,19)
    xticks(arange(16,20),['16','17','18','19'])
    figtext(0.6,0.4,"1E~2259+589",size=18)
    xlabel(r"$\lambda$ (\AA)")
    ylabel(r"$F_\lambda$ (ph s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)")
    
def plot_ne():
    os.chdir('/data/xmm/0142')
    l,f,e=xspecs.data_in('0142_1.fits')
    close()
    figure(1,(4.5,8.5))
    subplot(411)
    edge(l,f,e,"Ne",13,16)
    ylim(4e-4,2e-3)
    xticks(arange(13,17))
    os.chdir('../1048')
    l,f,e=xspecs.data_in('1048_1.fits')
    subplot(412)
    edge(l,f,e,"Ne",13,16)
    ylim(1e-5,5e-5)
    xticks(arange(13,17))
    ylabel(r"$F_\lambda$ (ph s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)")
    os.chdir('../2259')
    l,f,e=xspecs.data_in('2259_1.fits')
    subplot(413)
    edge(l,f,e,"Ne",13,16)
    ylim(9e-5,4.5e-4)
    xticks(arange(13,17))
    os.chdir('../1708')
    l,f,e=xspecs.data_in('1708_1.fits')
    subplot(414)
    edge(l,f,e,"Ne",13,16)
    ylim(1e-5,5e-5)
    xticks(arange(13,17),['13','14','15','16'])
    xlabel(r"$\lambda$ (\AA)")
    figtext(0.6,0.9,"4U~0142+61",size=16)
    figtext(0.6,0.71,"1E~1048.1$-$5937",size=16)
    figtext(0.6,0.47,"1E~2259+589",size=16)
    figtext(0.55,0.27,"1RXS J170849.0$-$400910",size=12)

def plot_si():
    os.chdir('/data/xmm/0142')
    l,f,e=xspecs.data_in('0142_1.fits')
    close()
    figure(1,(4.5,8.5))
    subplot(411)
    edge(l,f,e,"Si",6.2,7.5)
    xticks(arange(6.5,8,0.5))
    ylim(4.5e-3,9e-3)
    os.chdir('../1048')
    l,f,e=xspecs.data_in('1048_1.fits')
    subplot(412)
    edge(l,f,e,"Si",6.2,7.5)
    xticks(arange(6.5,8,0.5))
    ylim(2.6e-4,5.2e-4)
    ylabel(r"$F_\lambda$ (ph s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)")
    os.chdir('../2259')
    l,f,e=xspecs.data_in('2259_1.fits')
    subplot(413)
    edge(l,f,e,"Si",6.2,7.5)
    xticks(arange(6.5,8,0.5))
    ylim(1.5e-3,3e-3)
    os.chdir('../1708')
    l,f,e=xspecs.data_in('1708_1.fits')
    subplot(414)
    edge(l,f,e,"Si",6.2,7.5)
    xticks(arange(6.5,8,0.5),["6.5",'7.0','7.5'])
    ylim(9.5e-4,1.9e-3)
    xlabel(r"$\lambda$ (\AA)")
    figtext(0.6,0.9,"4U~0142+61",size=16)
    figtext(0.6,0.71,"1E~1048.1$-$5937",size=16)
    figtext(0.6,0.47,"1E~2259+589",size=16)
    figtext(0.55,0.27,"1RXS J170849.0$-$400910",size=12)

    
def plot_mg():
    os.chdir('/data/xmm/0142')
    l,f,e=xspecs.data_in('0142_1.fits')
    close()
    figure(1,(4.5,8.5))
    subplot(411)
    edge(l,f,e,"Mg",8.5,10.5)
    xticks(arange(8.5,11,0.5),[""])
    ylim(2.6e-3,1.3e-2)
    os.chdir('../1048')
    l,f,e=xspecs.data_in('1048_1.fits')
    subplot(412)
    edge(l,f,e,"Mg",8.5,10.5)
    xticks(arange(8.5,11,0.5),[""])
    ylim(.8e-4,4e-4)
    ylabel(r"$F_\lambda$ (ph s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)")
    os.chdir('../2259')
    l,f,e=xspecs.data_in('2259_1.fits')
    subplot(413)
    edge(l,f,e,"Mg",8.5,10.5)
    xticks(arange(8.5,11,0.5),[""])
    ylim(6e-4,3e-3)
    os.chdir('../1708')
    l,f,e=xspecs.data_in('1708_1.fits')
    subplot(414)
    edge(l,f,e,"Mg",8.5,10.5)
    xticks(arange(8.5,11,0.5),['8.5','9.0','9.5','10.0','10.5'])
    ylim(2e-4,1e-3)
    xlabel(r"$\lambda$ (\AA)")
    figtext(0.6,0.9,"4U~0142+61",size=16)
    figtext(0.6,0.71,"1E~1048.1$-$5937",size=16)
    figtext(0.6,0.47,"1E~2259+589",size=16)
    figtext(0.55,0.27,"1RXS J170849.0$-$400910",size=12)



