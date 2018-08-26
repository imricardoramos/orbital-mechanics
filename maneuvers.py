import numpy as np
from coordinates import kep2cart, kep2eq, eq2cart
from models import atmosDensity, cubesat
from scipy import integrate
from constants import constants

class Maneuvers:

  def __init__(self,com):
    self._PERTURBATION_ATMDRAG_ = 0
    self._PERTURBATION_J2_ = 0
    self._PERTURBATION_SOLARPRESS_ = 0

    self.current_r = 0 
    self.current_v = 0
    self.current_t = 0

    self.segments = []

    self.current_r, self.current_v = kep2cart(com)

    self.rTrace = self.current_r
    self.vTrace = self.current_v
    self.tTrace = self.current_t
    self.comTrace = np.array(com)
    self.equinoctialTrace = np.array(kep2eq(com))

  def addSegment(self,segment):
    self.segments.append(segment)

  def addPerturbation(self,perturbation):
    if perturbation == "atmosphere":
      self._PERTURBATION_ATMDRAG_ = 1
    if perturbation == "J2":
      self._PERTURBATION_J2_ = 1
    if perturbation == "solar_pressure":
      self._PERTURBATION_SOLARPRESS_ = 1

  def conwell(self, stateVector,t):
    r = stateVector[0:3]
    v = stateVector[3:6]
    z = np.linalg.norm(r)-constants.Re
    if t % (60*60*24) < 100:
      print("Day:"+str(t/60/60/24)+"\tHeight: "+str(z/1000)+" km")

    if(z > 100e3):
      p = 0
      if self._PERTURBATION_ATMDRAG_:
        #Atmospheric Drag
        vrel = v - np.cross(constants.wE,r)
        Fd = -0.5*atmosDensity(z/1000)*np.linalg.norm(vrel)*vrel*(1/cubesat.BC)
        p = p + Fd
      if self._PERTURBATION_J2_:
        J2 = 0.00108263

        zz = r[2]

        rNorm = np.linalg.norm(r)
        bigterm = np.array([1/rNorm*(5*zz**2/rNorm**2-1),
                            1/rNorm*(5*zz**2/rNorm**2-1),
                            1/rNorm*(5*zz**2/rNorm**2-3)])

        #J2 Gravity Gradient
        FJ2 = (3/2)*J2*constants.mu_E*constants.Re**2/rNorm**4*bigTerm*r
        p = p + FJ2

      if self._PERTURBATION_SOLARPRESS_:
        #Solar Radiation Pressure
        PSR = 4.56e-6
        #Absorbing Area
        As = np.pi*(cubesat.dimensions[0]/2)**2
        #Shadow function TODO
        nu = 1
        #Radiation Pressure
        FR = -nu*PSR*CR*As/cubesat.mass
        p = p + FR
      
      E = np.linalg.norm(v)**2/2-constants.mu_E/np.linalg.norm(r)
      if t % (60*60*24) < 100:
        print("Orbital Energy:"+str(E));
      #Differential Equations
      drdt = v
      dvdt = -constants.mu_E*(r/(np.linalg.norm(r)**3))+p

      return np.append(drdt,dvdt)
    else:
      return np.append(-r,-v)
  def betts(self,equinoctial_coordinates,t):

    p = equinoctial_coordinates[0]
    f = equinoctial_coordinates[1]
    g = equinoctial_coordinates[2]
    h = equinoctial_coordinates[3]
    k = equinoctial_coordinates[4]
    L = equinoctial_coordinates[5]

    q = 1+f*np.cos(L)+g*np.sin(L)
    s = np.sqrt(1+h**2+k**2)
    r,v = eq2cart(equinoctial_coordinates)
    z = np.linalg.norm(r) - constants.Re
    if t % (60*60*24) < 100:
      print("equinocitalElements:"+str(equinoctial_coordinates))
      print("Day:"+str(t/60/60/24)+"\tHeight: "+str(z/1000)+" km")

    #Perturbations in RSW Frame
    Delta = np.array([0,0,0])

    if self._PERTURBATION_ATMDRAG_:
      #Atmospheric Drag
      vrel = v - np.cross(constants.wE,r)
      Fd = np.array([0,-0.5*atmosDensity(z/1000)*np.linalg.norm(vrel)**2*(1/cubesat.BC),0])
      Delta = Delta + Fd


    A = np.array([[0, 2*p/q*np.sqrt(p/constants.mu_E),0],
                  [np.sqrt(p/constants.mu_E)*np.sin(L), np.sqrt(p/constants.mu_E)*1/q*((q+1)*np.cos(L)+f), -np.sqrt(p/constants.mu_E)*g/q*(h*np.sin(L)-k*np.cos(L))],
                  [-np.sqrt(p/constants.mu_E)*np.cos(L), np.sqrt(p/constants.mu_E)*(1/q)*((q+1)*np.sin(L)+g), np.sqrt(p/constants.mu_E)*(f/q)*(h*np.sin(L)-k*np.cos(L))],
                  [0, 0, np.sqrt(p/constants.mu_E)*s**2*np.cos(L)/(2*q)],
                  [0, 0, np.sqrt(p/constants.mu_E)*s**2*np.sin(L)/(2*q)],
                  [0, 0, np.sqrt(p/constants.mu_E)*(1/q)*(h*np.sin(L)-k*np.cos(L))]])
    b = np.array([0, 0, 0, 0, 0, np.sqrt(constants.mu_E*p)*(q/p)**2])

    dotx = np.matmul(A,Delta) + b
    return dotx

  def propagate(self,time):
    print("Propagating...from day ",self.current_t/60/60/24," to ",(self.current_t+time)/60/60/24)
    #Integrate
    y0 = np.append(self.current_r,self.current_v)
    t = np.linspace(self.current_t,self.current_t+time,time/60)
    y = integrate.odeint(self.conwell,y0,t)
    # Update traces
    self.tTrace = np.append(self.tTrace,t)
    self.rTrace = np.vstack((self.rTrace, y[:,0:3]))
    self.vTrace = np.vstack((self.vTrace, y[:,3:6]))
    # Update last values
    self.current_r = self.rTrace[-1,:]
    self.current_v = self.vTrace[-1,:]
    self.current_t = self.current_t+time

  def propagate2(self,time):
    #Integrate
    y0 = kep2eq(self.comTrace)
    t = np.linspace(self.current_t,self.current_t+time,time/60)
    y = integrate.odeint(self.betts,y0,t)
    # Update traces
    self.tTrace = np.append(self.tTrace,t)
    self.equinoctialTrace = np.vstack((self.equinoctialTrace,y))
    for idx,row in enumerate(self.equinoctialTrace):
      perc = idx/self.equinoctialTrace.shape[0]*100
      if perc-np.floor(perc) < 0.001:
        print("Perc:"+str(perc))
      r,v = eq2cart(row)
      self.rTrace = np.vstack((self.rTrace,r))
    self.current_t = self.current_t+time

  def impulsive_maneuver(self,dv):
    self.current_v = self.current_v+dv

