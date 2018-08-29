import numpy as np
class constants:
  # Universal Gravity
  G = 6.67e-11

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

  #SUN
  Ms = 1.989e30 
  mu_S = G*Ms
