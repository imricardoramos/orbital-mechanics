import numpy as np

# Universal Gravity
G = 6.67e-11
# Astronomical unit (m):
AU = 149597870691

#EARTH
# Earth Mass
Me = 5.97e24
#Earth Radius
Re = 6378e3
#Earth Angular Speed
wE = np.array([0,0,7.2921159e-5])
#Gravity Constants
mu_E = G*Me
#g0 constant (m/s)
g0 = 9.806

#MOON
Mm = 7.347e22 
mu_M = G*Mm
#Moon Radius
Rm = 1.737e3

#SUN
Ms = 1.989e30 
mu_S = G*Ms

#OTHER
#Universal Gas Constant R = 8.31432 J/K-mole
R = 8.31432
