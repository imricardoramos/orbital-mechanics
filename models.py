import numpy as np
from constants import constants
def atmosDensity(z):
    i=0
    h = [ 0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000];
    rho = [1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5, 1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9, 2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12, 1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15];
    H = [ 7.310, 6.427, 6.546, 7.360, 8.342, 7.583, 6.661, 5.927, 5.533, 5.703, 6.782, 9.973, 13.243, 16.322, 21.652, 27.974, 34.934, 43.342, 49.755, 54.513, 58.019, 60.980, 65.654, 76.377, 100.587, 147.203, 208.020];
    for j in range(0,27):
        if z >= h[j] and z < h[j+1]:
            i = j
    if z == 1000:
        i = 27
    density = rho[i]*np.exp(-(z-h[i])/H[i]);
    return density
def lunarPositionAlmanac2013(date):
  #Coefficients for Computing Lunar Position (Table 12.1 Curtis)
  a = [0, 6.29, -1.27, 0.66, 0.21, -0.19, -0.11]
  b = [218.32, 135.0, 259.3, 235.7, 269.9, 357.5, 106.5]
  c = [481267.881, 477198.87, -413335.36, 890534.22, 954397.74, 35999.05, 966404.03]
  d = [0, 5.13, 0.28, -0.28, -0.17, 0, 0]
  e = [0, 93.3, 220.2, 318.3, 217.6, 0, 0]
  f = [0, 483202.03, 960400.89, 6003.13, -407332.21, 0, 0]
  g = [0.9508, 0.0518, 0.0095, 0.0078, 0.0028, 0, 0]
  h = [0, 135.0, 259.3, 253.7, 269.9, 0, 0]
  k = [0, 477198.87, -413335.38, 890534.22, 954397.70, 0, 0]

  year = date.year
  month = date.month
  day = date.day
  UT = date.hour+date.minute/60
  #Julian Day
  J0 = 367*year-int(7*(year+int((month+9)/12))/4)+int(275*month/9)+day+1721013.5
  JD = J0+UT/24
  #Julian Centuries
  T0 = (JD-2451545)/36525
  #Obliquiy (in deg)
  epsilon = 23.439 - 0.0130042*T0
  epsilon = epsilon*np.pi/180
  #Lunar Ecliptic Longitude
  sums = 0
  for i in range(1,7):
    sums = sums + a[i]*np.sin(b[i]+c[i]*T0)
  lambd = b[0]+c[0]*T0+sums
  lambd = lambd*np.pi/180
  #Lunar Ecliptic Latitude
  delta = 0
  for i in range(1,5):
    delta = d[i]*np.sin(e[i]+f[i]*T0)
  delta = delta*np.pi/180
  #Horizontal Parallax
  sums = 0
  for i in range(1,5):
    sums = sums + g[i]*np.cos(h[i]+k[i]*T0)
  HP = g[0] + sums
  HP = HP*np.pi/180

  rm = constants.Re/np.sin(HP)
  rmVect = np.array([rm*np.cos(delta)*np.cos(lambd),
                     rm*np.cos(epsilon)*np.cos(delta)*np.sin(lambd) - np.sin(epsilon)*np.sin(delta),
                     rm*np.sin(epsilon)*np.cos(delta)*np.sin(lambd) + np.cos(epsilon)*np.sin(delta)])
  return rmVect
def solarPosition(date):
  year = date.year
  month = date.month
  day = date.day
  UT = date.hour+date.minute/60
  #Julian Day
  J0 = 367*year-int(7*(year+int((month+9)/12))/4)+int(275*month/9)+day+1721013.5
  JD = J0+UT/24

  # Astronomical unit (km):
  AU = 149597870.691

  # Julian days since J2000:
  n = JD - 2451545
  # Julian centuries since J2000:
  cy = n/36525

  # Mean anomaly (deg):
  M = 357.528 + 0.9856003*n
  M = M % 360
  M = M*np.pi/180
  # Mean longitude (deg):
  L = 280.460 + 0.98564736*n
  L = L % 360
  L = L*np.pi/180
  # Apparent ecliptic longitude (deg):
  lamda = L + 1.915*np.sin(M) + 0.020*np.sin(2*M)
  lamda = lamda % 360
  lamda = lamda*np.pi/180
  # Obliquity of the ecliptic (deg):
  eps = 23.439 - 0.0000004*n
  eps = eps*np.pi/180
  # Unit vector from earth to sun:
  u = np.array([np.cos(lamda),
                np.sin(lamda)*np.cos(eps),
                np.sin(lamda)*np.sin(eps)])
  # Distance from earth to sun (km):
  rS = (1.00014 - 0.01671*np.cos(M) - 0.000140*np.cos(2*M))*AU
  # Distance from earth to sun (m):
  rS = rS*1e3
  # Geocentric position vector (km):
  r_S = rS*u;
  
  return r_S

class Thruster:
  #ISP = 720s
  isp = 720
  #Trust = 1N
  thrust = 5e-2

  #Calculated mass flow rate (kg/s)
  massFlowRate = thrust/(isp*constants.g0)

class Spacecraft:
  def __init__(self,wetMass,dryMass,area):
    self.wetMass = wetMass
    self.dryMass = dryMass
    self.area = area
    self.thruster = Thruster()

  def BC(self,wetMass):
    Cd = 2.2
    return wetMass/(Cd*self.area)

class CurtisSat(Spacecraft):
  def __init__(self):
    m = 100
    A = np.pi/4
    BC = m/(Cd*A)

    Spacecraft.__init__(self,m,m,A)

class Cubesat(Spacecraft):
  def __init__(self):
    #Spacecraft Data
    totalMass = 3
    propellantMass = 1
    dimensions = [10e-2,10e-2,30e-2]
    A = dimensions[0]*dimensions[1]
    Spacecraft.__init__(self,totalMass,totalMass-propellantMass,A)
