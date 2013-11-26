emax = [.100,.284,.400,.532,.707,.867,1.303,
       1.840,2.471,3.210,4.038,7.111,8.331,10.0]
c0 = [17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,
        202.7,342.7,352.2,433.9,629.0,701.2]
c1 = [608.1,267.9,18.8,66.8,145.8,-380.6,169.3, 
         146.8,104.7,18.7,18.7,-2.4,30.9,25.2]
c2 = [-2150.,-476.1,4.3,-51.4,-61.1,294.0,-47.7,
         -31.5,-17.0,0.0,0.0,0.75,0.0,0.0]

def sigism(energy):
    """Interstellar absorption cross-section per hydrogen atom
    for solar abundances from Balucinska-Church & McGammon (1992)
    energy is in eV"""
    E=energy/1.e3    
#          (convert to keV)
    for i in range(1,14):
        if E < emax[i]: break
    result = (c0[i] + c1[i]*E + c2[i]*E**2)/E**3 * 1e-24
    return result

""" original code in Fortran:
      do 100 i=1,14  
        if (E .lt. Emax(i)) goto 200
100   continue  
      i=14
200   sigism=(c0(i)+c1(i)*E+c2(i)*E*E)/E**3 * 1.E-24
      return
      end"""
