import numpy as np
import constants

########### ATMOSPHERIC MODELS ############
def USSA76(z):
  """
  US Standard Atmosphere 1976. Exponential interpolation.

  Ref: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539.pdf
       https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/52703/versions/1/previews/atmosphere.m/index.html
       Orbital Mechanics For Engineers, Curtis 2013
  """
  z = z/1e3

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
def MSISE90(z,solarEpoch="low"):
  """
  Mass Spectrometer - Incoherent Scatter atmospheric model. Tabulated, low, mean and extremely high solar activity cases, linear interpolation.

  Ref: http://www.braeunig.us/space/atmos.htm
  """
  z = z/1e3

  altitudes = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900]
  if solarEpoch == 'low':
    densities = [1.17E+00, 9.48E-02, 4.07E-03, 3.31E-04, 1.69E-05, 5.77E-07, 1.70E-08, 2.96E-09, 9.65E-10, 3.90E-10, 1.75E-10, 8.47E-11, 4.31E-11, 2.30E-11, 1.27E-11, 7.22E-12, 4.21E-12, 2.50E-12, 1.51E-12, 9.20E-13, 5.68E-13, 3.54E-13, 2.23E-13, 1.42E-13, 9.20E-14, 6.03E-14, 4.03E-14, 2.75E-14, 1.93E-14, 1.39E-14, 1.03E-14, 7.90E-15, 6.24E-15, 5.06E-15, 4.21E-15, 3.58E-15, 3.09E-15, 2.70E-15, 2.39E-15, 2.13E-15, 1.91E-15, 1.73E-15, 1.56E-15, 1.42E-15, 1.30E-15, 1.18E-15]
  if solarEpoch == 'mean':
    densities = [1.17E+00, 9.49E-02, 4.07E-03, 3.31E-04, 1.68E-05, 5.08E-07, 1.80E-08, 3.26E-09, 1.18E-09, 5.51E-10, 2.91E-10, 1.66E-10, 9.91E-11, 6.16E-11, 3.94E-11, 2.58E-11, 1.72E-11, 1.16E-11, 7.99E-12, 5.55E-12, 3.89E-12, 2.75E-12, 1.96E-12, 1.40E-12, 1.01E-12, 7.30E-13, 5.31E-13, 3.88E-13, 2.85E-13, 2.11E-13, 1.56E-13, 1.17E-13, 8.79E-14, 6.65E-14, 5.08E-14, 3.91E-14, 3.04E-14, 2.39E-14, 1.90E-14, 1.53E-14, 1.25E-14, 1.03E-14, 8.64E-15, 7.32E-15, 6.28E-15, 5.46E-15]
  if solarEpoch == 'high':
    densities = [1.16E+00, 9.41E-02, 4.04E-03, 3.28E-04, 1.68E-05, 2.78E-07, 2.34E-08, 4.93E-09, 2.23E-09, 1.28E-09, 8.28E-10, 5.69E-10, 4.08E-10, 3.00E-10, 2.25E-10, 1.71E-10, 1.32E-10, 1.03E-10, 8.05E-11, 6.35E-11, 5.04E-11, 4.02E-11, 3.23E-11, 2.60E-11, 2.10E-11, 1.70E-11, 1.38E-11, 1.13E-11, 9.21E-12, 7.55E-12, 6.20E-12, 5.10E-12, 4.20E-12, 3.47E-12, 2.88E-12, 2.38E-12, 1.98E-12, 1.65E-12, 1.37E-12, 1.15E-12, 9.59E-13, 8.04E-13, 6.74E-13, 5.67E-13, 4.77E-13, 4.03E-13]
  for i,altitude in enumerate(altitudes):
    if altitude >= z:
      #Linear Interpolation
      k = (z-altitudes[i-1])/(altitudes[i]-altitudes[i-1])
      density = k*(densities[i]-densities[i-1])+densities[i-1]
      return density;

def HarrisPriester(z):
  """
  Harris-Priester atmospheric model. Tabulated, linear interpolation.

  Ref: David A Vallado - Fundamentals of Astrodynamics and Applications
  """
  z = z/1e3

  altitudes = [100, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 280, 290, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 760, 780, 800, 840, 880, 920, 960, 1000]
  minDensities = [4.974e-7, 2.490e-8, 8.377e-9, 3.899e-9, 2.122e-9, 1.163e-9, 8.008e-10, 5.283e-10, 3.617e-10, 2.557e-10, 1.839e-10, 1.341e-10, 9.949e-11, 7.488e-11, 5.709e-11, 4.403e-11, 2.697e-11, 2.139e-11, 1.708e-11, 1.099e-11, 7.214e-12, 4.824e-12, 3.274e-12, 2.249e-12, 1.558e-12, 1.091e-12, 7.701e-13, 5.474e-13, 3.916e-13, 2.819e-13, 2.042e-13, 1.488e-13, 1.092e-13, 8.070e-14, 6.012e-14, 4.519e-14, 3.430e-14, 2.620e-14, 2.043e-14, 1.607e-14, 1.036e-14, 8.496e-15, 7.069e-15, 4.680e-15, 3.200e-15, 2.210e-15, 1.560e-15, 1.150e-15]
  maxDensities = [4.974e-7, 2.490e-8, 8.710e-9, 4.059e-9, 2.215e-9, 1.344e-9, 8.758e-10, 6.010e-10, 4.297e-10, 3.162e-10, 2.396e-10, 1.853e-10, 1.455e-10, 1.157e-10, 9.308e-11, 7.555e-11, 5.095e-11, 4.226e-11, 3.526e-11, 2.511e-11, 1.819e-11, 1.337e-11, 9.955e-12, 7.492e-12, 5.684e-12, 4.355e-12, 3.362e-12, 2.612e-12, 2.042e-12, 1.605e-12, 1.267e-12, 1.005e-12, 7.997e-13, 6.390e-13, 5.123e-13, 4.131e-13, 3.325e-13, 2.691e-13, 3.325e-13, 1.779e-13, 1.190e-13, 9.776e-14, 8.059e-14, 5.741e-14, 4.210e-14, 3.130e-14, 2.360e-14, 1.810e-14] 
  for i,altitude in enumerate(altitudes):
    if altitude >= z:
      #Linear Interpolation
      k = (z-altitudes[i-1])/(altitudes[i]-altitudes[i-1])
      maxDensity = k*(maxDensities[i]-maxDensities[i-1])+maxDensities[i-1]
      minDensity = k*(minDensities[i]-minDensities[i-1])+minDensities[i-1]
      return minDensity, maxDensity;

def JacchiaRoberts(z):
  """
  ---INCOMPLETE---
  Jacchia-Roberts atmospheric model.
  Ref: http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?db_key=AST&bibcode=1971CeMec...4..368R&letter=0&classic=YES&defaultprint=YES&whole_paper=YES&page=368&epage=368&send=Send+PDF&filetype=.pdf
       https://github.com/JuliaSpace/SatelliteToolbox.jl/tree/2a6cbcb1cbf748ec9169dc851cb0a8b9b689f81b/src/earth/atmospheric_models/jr1971
       http://www.dem.inpe.br/~val/atmod/roberc.zip
  """
  #--Temperature Profile for Altitudes between 90 and 125 km--
  a = 444.3807
  b = 0.02385
  c = -392.8292
  d = -0.0021357
  # Exospheric Temperature (Jacchia 1970)
  # Tinf = ?
  # T(Z0) = 183 K
  T0 = 183
  # Temperature at inflection point
  Tx = a + b*Tinf + c*np.exp(d*Tinf)

  F = (Tx - T0)/35**4
  C = [-89284375.0, 3542400.0, -52687.5, 340.5, -0.8]
  cumSum = 0
  for i in range(0,5):
    cumSum += C[i]*Z**i

  T = Tx + F*cumSum
  #--Molecular Mass Profile for Altitudes between 90 and 125 km--
  a = np.array([28.15204, -8.5586e-2, 1.2840e-4, -1.0056e-5, -1.0210e-5, 1.5044e-6, 9.9826e-8])
  n = np.arange(0,6,1)
  M = a*(z-100)**n

  #--Atmospheric Density for Altitudes between 90 and 105 km--
  # M at 90 km
  M0 = a*(90-100)**n
  # Density(Z0 = 90 km) = 3.46e-9 gm/cm**3
  density0 = 3.46e-9
  # T(Z0 = 90 km) = 183 K
  T0 = 183
  # Parenthesis problem in k (See Ref)
  k = -constants.g0/constants.R*(Tx-T0)

  density = density0*(M*T0/M0*T)*F1**k*np.exp(k*F2)


  #--Atmospheric Density for Altitudes above 105 km--
  # Ref Table 1
  Q = [0.78110, 0.00934, 0.00001289,0.20995, 0.0, 0.0]
  M = [28.0134, 39.948, 4.0026, 31.9988, 15.9994, 1.00797]
  alpha = [0.0,0.0,-0.38,0.0,0.0,0.0]

  # --Atmospheric Density for Altitudes Above 125 km--
  d = np.zeros((6,))
  density = M*d

def mixAtmosphericModels(date,rhoLow,rhoHigh):
  decimalYear = date.year+date.month/12+date.day/12/30
  # Solar Activity peak every 11 years
  T = 11
  # 1 means max activity, 0 means lowest activity
  factor = (np.cos(2*np.pi*decimalYear/T-np.pi/6)+1)**2/4
  #factor = (1-np.cos(2*np.pi*decimalYear/T))/2
  rho = rhoHigh*factor+rhoLow*(1-factor)
  return rho

########### CELESTIAL MOTION ############
def lunarPositionAlmanac2013(date):
  #Coefficients for Computing Lunar Position (Table 12.1 Curtis)
  a = [0, 6.29, -1.27, 0.66, 0.21, -0.19, -0.11]
  b = [218.32, 135.0, 259.3, 235.7, 269.9, 357.5, 186.5]
  c = [481267.881, 477198.87, -413335.36, 890534.22, 954397.74, 35999.05, 966404.03]
  d = [0, 5.13, 0.28, -0.28, -0.17, 0, 0]
  e = [0, 93.3, 228.2, 318.3, 217.6, 0, 0]
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
    sums = sums + a[i]*np.sin((b[i]+c[i]*T0)*np.pi/180)
  lambd = b[0]+c[0]*T0+sums
  lambd = lambd % 360
  lambd = lambd*np.pi/180
  #Lunar Ecliptic Latitude
  delta = 0
  for i in range(1,5):
    delta = delta + d[i]*np.sin((e[i]+f[i]*T0)*np.pi/180)
  delta = delta % 360
  delta = delta*np.pi/180
  #Horizontal Parallax
  sums = 0
  for i in range(1,5):
    sums = sums + g[i]*np.cos((h[i]+k[i]*T0)*np.pi/180)
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

########### THRUSTERS ############
class Thruster:
  '''
  Generic Thruster Model
  '''
  def __init__(self,**kwargs):
    # DEFAULT VALUES
    self.name = "Generic Thruster"
    # Power [W]
    self.power = 1
    # ISP [s]
    self.isp = 720
    # Efficiency 
    self.eta = 1
    # Thrust [N]
    self.thrust = 1e-1

    for key in kwargs:
      if(key == 'thrust'):
        self.thrust = kwargs[key]
      if(key == 'isp'):
        self.isp = kwargs[key]
      if(key == 'power'):
        self.power = kwargs[key]
      if(key == 'throttle'):
        self.throttle = kwargs[key]
      if(key == 'efficiency'):
        self.eta = kwargs[key]

    self.massFlowRate = abs(self.thrust)/(self.isp*constants.g0)

  def operationalParams(self, Pa, throttle=1):
    # self.power is treated as max power the thruster can consume.
    # Pa = power available for the thruster to use
    # Pin = actual power consumed.
    Pa = max(Pa,0)
    Pin = min(Pa,throttle*self.power)
    #Power to thrust model
    thrust = self.thrust*Pin/self.power*self.eta
    massFlowRate = abs(thrust)/(self.isp*constants.g0)

    return Pin, thrust, massFlowRate

class IFMNanoThruster(Thruster):
  # Ref: https://www.cubesatshop.com/wp-content/uploads/2017/04/ENP-IFM-Nano-Thruster-Product-Overview.pdf 
  def __init__(self):
    Thruster.__init__(self, thrust=350e-6,isp=3000,throttle=0.5,efficiency=0.3)
    self.name = "IFM NanoTruster"

class NanoPropCGP3(Thruster):
  # Ref: https://gomspace.com/UserFiles/Subsystems/flyer/gomspace_nanoprop_cgp3.pdf
  # - Thrust: 1mN
  # - Thrust resolution: 10 uN
  # - Specific impulse: 60-110 sec 
  # - Total impulse 40 Ns
  # - Power consumption < 2W (average)
  # - Mass 300/350 gram (dry/wet)
  # - Operating pressure: 2-5 bar
  # - Temperature range  0° to 50°C 
  def __init__(self):
    Thruster.__init__(self, thrust=1e-3, isp=100, power=1.5)
    self.name = "NanoProp CGP3"

class VikiThruster(Thruster):
  def __init__(self):
    Thruster.__init__(self, thrust=4.565e-3, isp=6.8249, power=1.1794)
    self.name = "VikiThruster"

########### SOLAR PANELS ############
class solarPanels:
  '''
  Generic Solar Panels Model
  '''
  def __init__(self,n):
    self.name = "Generic Solar Panel"
    self.n = n
    # 10cm x 10cm
    self.area = 10e-2*10e-2
    self.efficiency = 0.2 
    #Power Density at Earth Surface (1.4 kW/m^2)
    self.nominalPowerDensity = 1400
    self.nominalPower = self.area*self.efficiency*self.nominalPowerDensity
  
  def power(self,r,datetime):
    uhat, r_s = solarPosition(datetime)
    r_sunSat = r-r_s

    #Shadow function
    rSNorm = np.linalg.norm(r_s)
    rNorm = np.linalg.norm(r)
    theta = np.arccos(np.dot(r_s,r)/(rSNorm*rNorm))
    theta1 = np.arccos(constants.Re/rNorm)
    theta2 = np.arccos(constants.Re/rSNorm)

    # SIGMOID CORRECTION
    if abs(theta1+theta2-theta) < 1*np.pi/180:
      nu = (np.sin((theta1+theta2-theta)/(1*np.pi/180)*np.pi/2)+1)/2
    else:
      if(theta1+theta2 > theta):
        nu = 1
      else:
        nu = 0

    outputPower = self.nominalPower*constants.AU**2/np.linalg.norm(r_sunSat)**2*nu
    return outputPower

class DHV_CS_10(solarPanels):
  # Ref: https://www.cubesatshop.com/product/cubesat-solar-panels/
  def __init__(self,n):
    self.name = "DHV CS 10"
    # Nominal Conditions: AM0 WRC = 1367 W/m2; T = 28 °C
    self.n = n
    self.area = 82.5e-3 * 98e-3
    self.efficiency = 0.3
    self.nominalPowerDensity = 1367
    self.nominalPower = 4.82*0.5*n

########### BATTERIES ############
class Battery:
  def __init__(self,P,S,cellVoltage,cellCapacity,chargeCurrent,dischargeCurrent):
    self.name = "Generic Battery"
    # Number of cells in parallel and series
    self.P = P
    self.S = S
    # Battery Voltage
    self.voltage = cellVoltage*S
    # Total capacity (Assuming cells are in parallel)
    self.capacity = cellCapacity*P
    # Battery Energy capacity [Wh]
    self.energy = self.voltage*self.capacity/1e3
    # Discharge Power
    self.dischargePower = self.voltage*dischargeCurrent*P/1e3
    # Charge Power
    self.chargePower = self.voltage*chargeCurrent*P/1e3

class NanoPowerBP4(Battery):
  # Ref: https://gomspace.com/UserFiles/Subsystems/datasheet/gs-ds-nanopower-bp4-27.pdf
  #      https://gomspace.com/UserFiles/Subsystems/datasheet/gs-ds-nanopower-battery-17.pdf
  def __init__(self,configuration):
    if configuration == "2P-2S":
      Battery.__init__(self, 2, 2, 3.7, 2600, 1000, 1000)
    if configuration == "1P-4S":
      Battery.__init__(self, 1, 4, 3.7, 2600, 1000, 1000)
    self.name = "NanoPower BP4"

#class BA0X(Battery):
#   Ref: https://www.cubesatshop.com/wp-content/uploads/2016/11/EXA-BA0x-Brochure.pdf
#   <!>NO CHARGE/DISCHARGE CURRENT INFO<!>
#
#  def __init__(self,model):
#    if model == "BA01/S":
#      Battery.__init__(self,8,3.7,900)
#    if model == "BA01/D":
#      Battery.__init__(self,16,3.7,900)
#    if model == "BA02/S":
#      Battery.__init__(self,6,3.7,900)
#    if model == "BA02/D":
#      Battery.__init__(self,12,3.7,900)
#    self.name = model

########### SPACECRAFTS ############
class Spacecraft:
  def __init__(self,wetMass,dryMass,area):
    self.wetMass = wetMass
    self.dryMass = dryMass
    self.area = area
    self.Cd = 2.2
    self.Cr = 2

    self.thruster = Thruster()
    self.solarPanels = solarPanels(1)
    self.battery = Battery(8,1,3.7,900,1000,1000)

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

