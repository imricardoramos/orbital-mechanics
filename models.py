import numpy as np
import constants
#def getScaleHeight(h):
#    kms = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550];
#    scaleHeights = [8.4, 5.9, 25.5, 37.5, 44.8, 50.3, 54.8, 58.2, 61.3, 64.5, 68.7];
#    for i,km in enumerate(kms):
#        if h < km:
#            #interpolation
#            k = (kms[i]-h)/(kms[i]-kms[i-1])
#            scaleHeight = k*(scaleHeights[i]-scaleHeights[i-1])+scaleHeights[i-1]
#            return scaleHeight;
#def getAtmDensity(h,solarEpoch='low'):
#    kms = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520];
#    if solarEpoch == 'high':
#        densities = [1.16E+00, 9.41E-02, 4.04E-03, 3.28E-04, 1.68E-05, 2.78E-07, 2.34E-08, 4.93E-09, 2.23E-09, 1.28E-09, 8.28E-10, 5.69E-10, 4.08E-10, 3.00E-10, 2.25E-10, 1.71E-10, 1.32E-10, 1.03E-10, 8.05E-11, 6.35E-11, 5.04E-11, 4.02E-11, 3.23E-11, 2.60E-11, 2.10E-11, 1.70E-11, 1.38E-11];
#    elif solarEpoch == 'mid':
#        densities = [1.17E+00, 9.49E-02, 4.07E-03, 3.31E-04, 1.68E-05, 5.08E-07, 1.80E-08, 3.26E-09, 1.18E-09, 5.51E-10, 2.91E-10, 1.66E-10, 9.91E-11, 6.16E-11, 3.94E-11, 2.58E-11, 1.72E-11, 1.16E-11, 7.99E-12, 5.55E-12, 3.89E-12, 2.75E-12, 1.96E-12, 1.40E-12, 1.01E-12, 7.30E-13, 5.31E-13];
#    elif solarEpoch == 'low':
#        densities = [1.17E+00, 9.48E-02, 4.07E-03, 3.31E-04, 1.69E-05, 5.77E-07, 1.70E-08, 2.96E-09, 9.65E-10, 3.90E-10, 1.75E-10, 8.47E-11, 4.31E-11, 2.30E-11, 1.27E-11, 7.22E-12, 4.21E-12, 2.50E-12, 1.51E-12, 9.20E-13, 5.68E-13, 3.54E-13, 2.23E-13, 1.42E-13, 9.20E-14, 6.03E-14, 4.03E-14];
#    for i,km in enumerate(kms):
#        if h <= km:
#            #interpolation
#            k = (kms[i]-h)/(kms[i]-kms[i-1])
#            density = k*(densities[i]-densities[i-1])+densities[i-1]
#            return density;

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
  f = [0, 483202.03, 960400.89, 6003.15, -407332.21, 0, 0]
  g = [0.9508, 0.0518, 0.0095, 0.0078, 0.0028, 0, 0]
  h = [0, 135.0, 259.3, 253.7, 269.9, 0, 0]
  k = [0, 477198.87, -413335.38, 890534.22, 954397.70, 0, 0]

  year = date.year
  month = date.month
  day = date.day
  UT = date.hour+date.minute/60+date.second/60/60
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
  lambd = lambd % 360
  lambd = lambd*np.pi/180
  #Lunar Ecliptic Latitude
  delta = 0
  for i in range(1,5):
    delta = delta + d[i]*np.sin(e[i]+f[i]*T0)
  delta = delta % 360
  delta = delta*np.pi/180
  #Horizontal Parallax
  sums = 0
  for i in range(1,5):
    sums = sums + g[i]*np.cos(h[i]+k[i]*T0)
  HP = g[0] + sums
  HP = HP % 360
  HP = HP*np.pi/180

  rm = constants.Re/np.sin(HP)
  rmVect = rm*np.array([np.cos(delta)*np.cos(lambd),
                        np.cos(epsilon)*np.cos(delta)*np.sin(lambd) - np.sin(epsilon)*np.sin(delta),
                        np.sin(epsilon)*np.cos(delta)*np.sin(lambd) + np.cos(epsilon)*np.sin(delta)])
  return rmVect

def solarPosition(date,unit="meters"):
  year = date.year
  month = date.month
  day = date.day
  UT = date.hour+date.minute/60
  # Astronomical unit (km):
  AU = 149597870.691
  #Julian Day
  J0 = 367*year-int(7*(year+int((month+9)/12))/4)+int(275*month/9)+day+1721013.5
  JD = J0+UT/24

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
  # Apparent ecliptic longitude (deg):
  lamda = L + 1.915*np.sin(M) + 0.020*np.sin(2*M)
  L = L*np.pi/180
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
  rS = (1.00014 - 0.01671*np.cos(M) - 0.000140*np.cos(2*M))
  
  if unit == "AU":
    scaler = 1
  elif unit == "km":
    scaler = AU
  else:
    scaler = AU*1e3

  # Distance from earth to sun:
  rS = rS*scaler
  # Geocentric position vector:
  r_S = rS*u;
  
  return u, r_S

class Thruster:
  def __init__(self,**kwargs):
    # DEFAULT VALUES
    # Power Available [W]
    self.Pa = 1
    # ISP [s]
    self.isp = 720
    # Efficiency 
    self.eta = 1
    # Duty Cycle
    self.D = 1
    # Thrust [N]
    self.thrust = 1e-1

    for key in kwargs:
      if(key == 'thrust'):
        self.thrust = kwargs[key]
      if(key == 'isp'):
        self.isp = kwargs[key]
      if(key == 'power_available'):
        self.Pa = kwargs[key]
      if(key == 'duty_cycle'):
        self.D = kwargs[key]
      if(key == 'efficiency'):
        self.eta = kwargs[key]

    self.massFlowRate = abs(self.thrust)/(self.isp*constants.g0)

  def powerToThrust(self):
    #Power to thrust model
    self.thrust = 2*self.D*self.eta*self.Pa/(constants.g0*self.isp)
  def thrustToPower(self):
    return self.thrust*(constants.g0*self.isp)/(2*self.D*self.eta)

class IFMNanoThruster(Thruster):
  # Ref: https://www.cubesatshop.com/wp-content/uploads/2017/04/ENP-IFM-Nano-Thruster-Product-Overview.pdf 
  def __init__(self):
    Thruster.__init__(self, thrust=350e-6,isp=3000,duty_cycle=0.5,efficiency=0.3)

class solarPanels:
  def __init__(self):
    self.area = 30e-2*10e-2
    self.efficiency = 1
    #Power Density at Earth Surface (1.4 kW/m^2)
    self.nominalPowerDensity = 1400
    self.nominalPower = self.area*self.efficiency*self.nominalPowerDensity
  
  def power(self,r,datetime):
    uhat, r_s = solarPosition(datetime)
    r_sunSat = r-r_s
    outputPower = self.nominalPower*constants.AU**2/np.linalg.norm(r_sunSat)**2
    return outputPower

class DHV_CS_10(solarPanels):
  # Ref: https://www.cubesatshop.com/product/cubesat-solar-panels/
  def __init__(self,n):
    # Nominal Conditions: AM0 WRC = 1367 W/m2; T = 28 °C
    self.area = 82.5e-3 * 98e-3 * 0.8
    self.efficiency = 0.3
    self.nominalPowerDensity = 1367
    self.nominalPower = 4.82*0.5*n

class Spacecraft:
  def __init__(self,wetMass,dryMass,area):
    self.wetMass = wetMass
    self.dryMass = dryMass
    self.area = area
    self.solarPanels = solarPanels()
    self.Cd = 2.2

  def BC(self,wetMass):
    return wetMass/(self.Cd*self.area)

class CurtisSat(Spacecraft):
  def __init__(self):
    m = 100
    A = np.pi/4
    Spacecraft.__init__(self,m,m,A)

class Cubesat(Spacecraft):
  def __init__(self,nUnits):
    if nUnits == "1U":
      totalMass = 1
      propellantMass = totalMass/3
      dimensions = [10e-2,10e-2,30e-2]
    if nUnits == "3U":
      totalMass = 3
      propellantMass = totalMass/3
      dimensions = [10e-2,10e-2,30e-2]
    if nUnits == "6U":
      totalMass = 6
      propellantMass = totalMass/3
      dimensions = [10e-2,10e-2,30e-2]
    if nUnits == "12U":
      totalMass = 12
      propellantMass = totalMass/3
      dimensions = [10e-2,10e-2,30e-2]

    A = 2*(dimensions[0]*dimensions[1]+\
           dimensions[1]*dimensions[2]+\
           dimensions[0]*dimensions[2])/6

    Spacecraft.__init__(self,totalMass,totalMass-propellantMass,A)

