# function from Allison Smith's github:
# https://github.com/kallisons/pO2_conversion/function_pO2.py
# with some modifications for style and units

# pO2 Conversion
# Oxygen concentration is converted to percent oxygen saturation using the equations from Garcia and Gordon (1992).
# The percent oxygen saturation is divided by 0.21 (the fractional atmospheric concentration of oxygen) to get pO2 in atmospheres (atm).
# pO2 is then corrected for the hydrostatic pressure at depth (Enns et al., 1965).
# The units for pO2 are converted to kilopascals (kPa), the SI Units for pressure.
# References:
# - García HE, Gordon LI (1992) Oxygen solubility in seawater: Better fitting equations. Limnology and Oceanography, 37, 1307–1312.
# - Enns T, Scholander PF, Bradstreet ED (1965) Effect of hydrostatic pressure on gases dissolved in water. The Journal of Physical Chemistry, 69, 389–391.
 
##UNITS
#o2 in umol/Kg
#temp in Celsius (= potential temperature, NOT in situ)
#sal in psu
#depth in m

import numpy as np

def calc_po2(o2, temp, sal, depth):
    """Computes po2 from o2 [umol/kg], potential temperature [Celsius], salinity [psu], depth [m]."""
    a_0 = 5.80871
    a_1 = 3.20291
    a_2 = 4.17887
    a_3 = 5.10006
    a_4 = -9.86643e-2
    a_5 = 3.80369
    b_0 = -7.01577e-3
    b_1 = -7.70028e-3
    b_2 =  -1.13864e-2
    b_3 = -9.51519e-3
    c_0 = -2.75915E-7

    tt = 298.15 - temp
    tk = 273.15 + temp
    ts = np.log(tt / tk)

    #correct for pressure at depth
    V = 32e-6 #partial molar volume of O2 (m3/mol)
    R = 8.31 #Gas constant [J/mol/K]
    db2Pa = 1e4 #convert pressure: decibar to Pascal
    atm2Pa = 1.01325e5 #convert pressure: atm to Pascal

    #calculate pressure in dB from depth in m
    pres = depth*(1.0076+depth*(2.3487e-6 - depth*1.2887e-11));

    #convert pressure from decibar to pascal
    dp = pres*db2Pa
    pCor = np.exp((V*dp)/(R*(temp+273.15)))

    o2_sat = np.exp(a_0 + a_1*ts + a_2*ts**2 + a_3*ts**3 + a_4*ts**4 + a_5*ts**5 + sal*(b_0 + b_1*ts + b_2*ts**2 + b_3*ts**3) + c_0*sal**2)
    
    o2_alpha = (o2_sat / 0.21)  #0.21 is atmospheric composition of O2
    kh = o2_alpha*pCor
    po2 = (o2 / kh)*101.32501  #convert po2 from atm to kPa

    return po2
