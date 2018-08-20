import numpy as np
from coordinates import kep2cart
from models import atmosDensity, cubesat
from scipy import integrate
from constants import constants

class Maneuvers:
  _PERTURBATION_ATMDRAG_ = 0
  _PERTURBATION_J2_ = 0
  _PERTURBATION_RADIATION_ = 0

  current_r = 0 
  current_v = 0
  current_t = 0

  rTrace = np.array([])
  vTrace = np.array([])
  tTrace = np.array([])

  com = []
  segments = []

  def __init__(self,com):
    self.com = com
    self.current_r, self.current_v = kep2cart(com)
    self.rTrace = self.current_r
    self.vTrace = self.current_v
    self.tTrace = self.current_t

  def addSegment(self,segment):
    self.segments.append(segment)

  def conwell(self, stateVector,t):
    r = stateVector[0:3]
    v = stateVector[3:6]
    z = np.linalg.norm(r)-constants.Re
    if t % (60*60*24) < 100:
      print("Day:"+str(t/60/60/24)+"\tHeight: "+str(z/1000)+" km")

    if(z > 100e3):
      vrel = v - np.cross(constants.wE,r)
      p = 0
      if _PERTURBATIONS_ATMDRAG_:
        #Atmospheric Drag
        Fd = -0.5*atmosDensity(z/1000)*np.linalg.norm(vrel)*vrel*(1/cubesat.BC)
        p = p + Fd
      if _PERTURBATIONS_J2_:
        J2 = 0.00108263
        #TODO
        phi=0

        rNorm = np.linalg.norm(r)
        z = rNorm*np.cos(phi)
        bigterm = np.array([1/rNorm*(5*z**2/rNorm**2-1),1/rNorm*(5*z**2/rNorm**2-1),1/rNorm*(5*z**2/rNorm**2-3)])

        #J2 Gravity Gradient
        FJ2 = (3/2)*J2*constants.mu_E*constants.Re**2/rNorm**4*bigTerm*r
        p = p + FJ2

      if _PERTURBATIONS_RAD_:
        #Solar Radiation Pressure
        PSR = 4.56e-6
        #Absorbing Area
        As = np.pi*(cubesat.dimensions[0]/2)**2
        #Shadow function TODO
        nu = 1
        #Radiation Pressure
        FR = -nu*PSR*CR*As/cubesat.mass
        p = p + FR
      
      #Differential Equations
      drdt = v
      dvdt = -constants.mu_E*r/(np.linalg.norm(r)**3)+p

      return np.append(drdt,dvdt)
    else:
      return np.append(-r,-v)

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

  def impulsive_maneuver(self,dv):
    self.current_v = self.current_v+dv

