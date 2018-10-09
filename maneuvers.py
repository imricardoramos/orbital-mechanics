import numpy as np
from coordinates import kep2cart, kep2eq, eq2cart, eq2kep
import models
from scipy import integrate
from constants import constants
from datetime import timedelta

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as  plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import mplcursors

import auxiliary

class ManeuversHistory:
  def __init__(self):
    self.t = [0]
    self.datetime = None 
    self.r = None
    self.v = None
    self.coe = None
    self.mee = None
    self.propMass = None
    self.dv = None
    self.maneuverIdxs = [0]

class Maneuvers:

  def __init__(self,coe,spacecraft,startDate,**kwargs):

    #Perturbation Flags
    self._PERTURBATION_ATMDRAG_ = 0
    self._PERTURBATION_J2_ = 0
    self._PERTURBATION_SOLARPRESS_ = 0
    self._PERTURBATION_THRUST_ = 0
    self._PERTURBATION_MOON_ = 0
    self._PERTURBATION_SUN_ = 0

    self.spacecraft = spacecraft
    self.startDate = startDate

    # Set up history to store data along maneuvers
    self.history = ManeuversHistory()
    self.history.r, self.history.v = kep2cart(coe)
    self.history.coe = np.array(coe)
    self.history.mee = kep2eq(coe)
    self.history.mee.shape = (1,6)
    self.history.propMass = [self.spacecraft.wetMass-self.spacecraft.dryMass]
    self.history.datetime = startDate

    # Other options
    self.verbose = False;
    for key in kwargs:
      if key == 'verbose':
        self.verbose = kwargs[key]

  def addPerturbation(self,perturbation):
    if perturbation == "atmosphere":
      self._PERTURBATION_ATMDRAG_ = 1
    if perturbation == "J2":
      self._PERTURBATION_J2_ = 1
    if perturbation == "solar_pressure":
      self._PERTURBATION_SOLARPRESS_ = 1
    if perturbation == "thrust":
      self._PERTURBATION_THRUST_ = 1
    if perturbation == "moon_gravity":
      self._PERTURBATION_MOON_ = 1
    if perturbation == "sun_gravity":
      self._PERTURBATION_SUN_ = 1

  def removePerturbation(self,perturbation):
    if perturbation == "atmosphere":
      self._PERTURBATION_ATMDRAG_ = 0
    if perturbation == "J2":
      self._PERTURBATION_J2_ = 0
    if perturbation == "solar_pressure":
      self._PERTURBATION_SOLARPRESS_ = 0
    if perturbation == "thrust":
      self._PERTURBATION_THRUST_ = 0
    if perturbation == "moon_gravity":
      self._PERTURBATION_MOON_ = 0
    if perturbation == "sun_gravity":
      self._PERTURBATION_SUN_ = 0

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
        Fd = -0.5*models.atmosDensity(z/1000)*np.linalg.norm(vrel)*vrel*(1/self.spacecraft.BC(self.spacecraft.dryMass+propMass))
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
        As = self.spacecraft.area
        #Shadow function TODO
        nu = 1
        #Radiation Pressure
        FR = -nu*PSR*CR*As/self.spacecraft.mass
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
  def betts(self,y,t):
    equinoctial_coordinates = y[0:6]
    p = equinoctial_coordinates[0]
    f = equinoctial_coordinates[1]
    g = equinoctial_coordinates[2]
    h = equinoctial_coordinates[3]
    k = equinoctial_coordinates[4]
    L = equinoctial_coordinates[5]

    propMass = y[6]

    q = 1+f*np.cos(L)+g*np.sin(L)
    s = np.sqrt(1+h**2+k**2)
    r,v = eq2cart(equinoctial_coordinates)
    z = np.linalg.norm(r) - constants.Re
    if self.verbose:
      if t % (60*60*24) < 100:
        #print("equinocitalElements:"+str(equinoctial_coordinates))
        print("Day:"+str(t/60/60/24)+"\tHeight: "+str(z/1000)+" km"+"\tMass: "+str(propMass))

    #RSW frame definition from r and v
    rhat = r/np.linalg.norm(r)
    w = np.cross(rhat,v)
    what = w/np.linalg.norm(w)
    shat = np.cross(what,rhat)
    rsw = [rhat, shat, what]

    #Perturbations in RSW Frame
    Delta = np.array([0,0,0])

    if self._PERTURBATION_ATMDRAG_:
      #Atmospheric Drag
      vrel = v - np.cross(constants.wE,r)
      Fd = -0.5*models.atmosDensity(z/1000)*np.linalg.norm(vrel)*vrel*(1/self.spacecraft.BC(self.spacecraft.dryMass+propMass))
      Delta = Delta + np.dot(rsw,Fd)

    if self._PERTURBATION_J2_:
      J2 = 0.00108263

      zz = r[2]

      rNorm = np.linalg.norm(r)
      bigTerm = np.array([1/rNorm*(5*zz**2/rNorm**2-1),
                          1/rNorm*(5*zz**2/rNorm**2-1),
                          1/rNorm*(5*zz**2/rNorm**2-3)])

      #J2 Gravity Gradient
      FJ2 = (3/2)*J2*constants.mu_E*constants.Re**2/rNorm**4*bigTerm*r
      Delta = Delta + np.dot(rsw,FJ2)

    if self._PERTURBATION_SOLARPRESS_:
      #u-hat vector pointing from Sun to Earth and Sun position vector
      uhat, rS = models.solarPosition(self.startDate+timedelta(seconds=t))
      
      #Solar Radiation Pressure
      PSR = 4.56e-6
      #Absorbing Area
      As = self.spacecraft.area

      #Shadow function
      rSNorm = np.linalg.norm(rS)
      rNorm = np.linalg.norm(r)
      theta = np.arccos(np.dot(rS,r)/(rSNorm*rNorm))
      theta1 = np.arccos(constants.Re/rNorm)
      theta2 = np.arccos(constants.Re/rSNorm)
      if(theta1+theta2 > theta):
        nu = 1
      else:
        nu = 0

      #Radiation Pressure Coefficient (lies between 0 and 1)
      CR = 0.5
      #Spacecraft mass
      mass = self.spacecraft.dryMass+propMass
      #Radiation Pressure acceleration
      FR = -nu*PSR*CR*As/mass*uhat
      Delta = Delta + np.dot(rsw,FR)
      #print(Delta)

    if self._PERTURBATION_MOON_:
      r_m = models.lunarPositionAlmanac2013(self.startDate+timedelta(minutes=t))
      r_ms = r_m-r 
      pMoon = constants.mu_M*(r_ms/np.linalg.norm(r_ms)**3-r_m/np.linalg.norm(r_m)**3)
      Delta = Delta + np.dot(rsw,pMoon)

    if self._PERTURBATION_SUN_:
      uhat, r_s = models.solarPosition(self.startDate+timedelta(seconds=t))
      r_sunSat = r_s-r 
      #F(q) formula from F.3 Appendix Curtis 2013
      # 	c = b - a; a << b
      # 	F = 1 - c**3/b**2
      qq = np.dot(r,(2*r_s-r))/np.linalg.norm(r_s)**2
      Fq = (qq**2 - 3*qq + 3)*qq/(1+(1-qq)**(3/2))
      pSun = constants.mu_S/np.linalg.norm(r_sunSat)**3*(Fq*r_s-r)
      Delta = Delta + np.dot(rsw,pSun)

    if self._PERTURBATION_THRUST_ and propMass > 0:
      mass = self.spacecraft.dryMass+propMass
      #Thrust in out-of-plane direction
      #Delta = Delta + [0,0,self.spacecraft.thruster.thrust/mass]
      #Thrust in speed direction
      Fth = self.spacecraft.thruster.thrust*v/np.linalg.norm(v)
      Delta = Delta + np.dot(rsw,Fth/mass)

      dpropMass = -self.spacecraft.thruster.massFlowRate
    else:
      dpropMass = 0


    A = np.array([[0, 2*p/q*np.sqrt(p/constants.mu_E),0],
                  [np.sqrt(p/constants.mu_E)*np.sin(L), np.sqrt(p/constants.mu_E)*1/q*((q+1)*np.cos(L)+f), -np.sqrt(p/constants.mu_E)*g/q*(h*np.sin(L)-k*np.cos(L))],
                  [-np.sqrt(p/constants.mu_E)*np.cos(L), np.sqrt(p/constants.mu_E)*(1/q)*((q+1)*np.sin(L)+g), np.sqrt(p/constants.mu_E)*(f/q)*(h*np.sin(L)-k*np.cos(L))],
                  [0, 0, np.sqrt(p/constants.mu_E)*s**2*np.cos(L)/(2*q)],
                  [0, 0, np.sqrt(p/constants.mu_E)*s**2*np.sin(L)/(2*q)],
                  [0, 0, np.sqrt(p/constants.mu_E)*(1/q)*(h*np.sin(L)-k*np.cos(L))]])
    b = np.array([0, 0, 0, 0, 0, np.sqrt(constants.mu_E*p)*(q/p)**2])

    dotx = np.matmul(A,Delta) + b
    return np.append(dotx,dpropMass)

  def propagate(self,time):
    print("Propagating...from day ",self.history.t[-1]/60/60/24," to ",(self.history.t[-1]+time)/60/60/24)
    #Integrate
    y0 = np.append(self.history.r[-1,:],self.history.v[-1,:])
    t = np.linspace(self.history.t[-1],self.history.t[-1]+time,time/60)
    y = integrate.odeint(self.conwell,y0,t)
    # Update history
    self.history.t = np.append(self.history.t,t)
    self.history.r = np.vstack((self.history.r, y[:,0:3]))
    self.history.v = np.vstack((self.history.v, y[:,3:6]))

  def propagate2(self,time):
    print("Propagating...from day ",self.history.t[-1]/60/60/24," to ",(self.history.t[-1]+time)/60/60/24)
    #Integrate
    y0 = np.append(self.history.mee[-1,:],self.history.propMass[-1])
    t = np.linspace(self.history.t[-1],self.history.t[-1]+time,time/60)
    sol = integrate.odeint(self.betts,y0,t)
    y = sol[:,0:6] 
    propMass = sol[:,6]

    #Initialize r and v for faster stacking
    rLocalHist = np.zeros([y.shape[0],3])
    vLocalHist = np.zeros([y.shape[0],3])
    #Initialize coe faster stacking
    coeLocalHist = np.zeros([y.shape[0],6])

    for idx,row in enumerate(y):
      perc = idx/y.shape[0]*100
      if self.verbose:
        if perc % 10 == 0:
          print(str(perc)+"%")
      r,v = eq2cart(row)
      rLocalHist[idx,:] = r
      vLocalHist[idx,:] = v
      coeLocalHist[idx,:] = eq2kep(row)
    datetime = np.array([self.startDate + timedelta(seconds=seconds) for seconds in t])

    # Update history
    self.history.t   = np.append(self.history.t,t)
    self.history.mee = np.vstack((self.history.mee,y))
    self.history.coe = np.vstack((self.history.coe,coeLocalHist))
    self.history.r   = np.vstack((self.history.r,rLocalHist))
    self.history.v   = np.vstack((self.history.v,vLocalHist))
    self.history.maneuverIdxs.append(len(self.history.t))
    self.history.propMass = np.append(self.history.propMass,propMass)
    self.history.datetime = np.append(self.history.datetime,datetime)

  def propagate3(self,time):
    print("Propagating...from day ",self.history.t[-1]/60/60/24," to ",(self.history.t[-1]+time)/60/60/24)
    #Integrate
    y0 = np.append(self.history.mee[-1,:],self.history.propMass[-1])
    t = np.linspace(self.history.t[-1],self.history.t[-1]+time,time/60)
    sol = integrate.solve_ivp(self.betts,(t[0],t[-1]),y0,method="RK45")
    t = sol.t
    y = sol.y[:,0:6] 
    propMass = sol.y[:,6]

    #Initialize r and v for faster stacking
    rLocalHist = np.zeros([y.shape[0],3])
    vLocalHist = np.zeros([y.shape[0],3])

    for idx,row in enumerate(y):
      perc = idx/y.shape[0]*100
      if perc % 10 == 0:
        print(str(perc)+"%")
      r,v = eq2cart(row)
      rLocalHist[idx,:] = r
      vLocalHist[idx,:] = v
    datetime = np.array([self.startDate + timedelta(seconds=seconds) for seconds in t])

    # Update history
    self.history.t   = np.append(self.history.t,t)
    self.history.mee = np.vstack((self.history.mee,y))
    self.history.r   = np.vstack((self.history.r,rLocalHist))
    self.history.v   = np.vstack((self.history.v,vLocalHist))
    self.history.maneuverIdxs.append(len(self.history.t))
    self.history.propMass = np.append(self.history.propMass,propMass)
    self.history.datetime = np.append(self.history.datetime,datetime)
  def impulsive_maneuver(self,dv):
    self.current_v = self.current_v+dv

  def plot(self, item, itemHistory=None):
    if item == "com":
      # PLOTTING CLASSICAL ORBITAL ELEMENTS
      titles = ["a","e","i","$\omega$","$\Omega$","M"]
      ylabels = ["[m]", "", "[°]", "[°]", "[°]", "[°]"]
      fig, axes = plt.subplots(3,2,figsize=(10,8))
      for i in range(0,6):
        for j in range(0,len(self.history.maneuverIdxs)-1):
          maneuverSlice = slice(self.history.maneuverIdxs[j],self.history.maneuverIdxs[j+1])
          if i in [2,3,4,5]:
            axes[int((i-i%2)/2),i%2].plot(self.history.datetime[maneuverSlice],self.history.coe[maneuverSlice,i]*180/np.pi)
          else:
            axes[int((i-i%2)/2),i%2].plot(self.history.datetime[maneuverSlice],self.history.coe[maneuverSlice,i])
          axes[int((i-i%2)/2),i%2].set_title(titles[i]+" "+ylabels[i])
          
          fig.autofmt_xdate()
          axes[int((i-i%2)/2),i%2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
          axes[int((i-i%2)/2),i%2].yaxis.get_major_formatter().set_scientific(False)
          axes[int((i-i%2)/2),i%2].yaxis.get_major_formatter().set_useOffset(False)
          axes[int((i-i%2)/2),i%2].grid(b=True)

    if item == "3d-trajectory":
      #Plot 3D Trajectory
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      markers = np.zeros([len(self.history.maneuverIdxs)-1,3])
      
      for i in range(0,len(self.history.maneuverIdxs)-1):
          maneuverSlice = slice(self.history.maneuverIdxs[i],self.history.maneuverIdxs[i+1])
          ax.plot3D(self.history.r[maneuverSlice,0]/1000,
                    self.history.r[maneuverSlice,1]/1000,
                    self.history.r[maneuverSlice,2]/1000,linewidth=1)
          markers[i,:]= self.history.r[self.history.maneuverIdxs[i],:]/1000
      ax.plot3D(markers[:,0],markers[:,1],markers[:,2],"k.")
      auxiliary.set_axes_equal(ax)
      ax.set_aspect("equal")
      scale_x = 1.2
      scale_y = 1.2
      scale_z = 1.2
      ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1]))
      ax.set_title("Satellite Trajectory [km]")
      ax.set_xlabel("X [km]");
      ax.set_ylabel("Y [km]");
      ax.set_zlabel("Z [km]");


    if item == "singleItem":
      if np.isscalar(itemHistory) and itemHistory == None:
        raise Exception("History Data not specified.")
      else:
        fig, ax = plt.subplots(figsize=(10,4))
        for i in range(0,len(self.history.maneuverIdxs)-1):
            maneuverSlice = slice(self.history.maneuverIdxs[i],self.history.maneuverIdxs[i+1])
            ax.plot(self.history.datetime[maneuverSlice], itemHistory[maneuverSlice],linewidth=1)

        ax.set_title("Satellite Height [km]")
        fig.autofmt_xdate()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.yaxis.get_major_formatter().set_scientific(False)
        ax.yaxis.get_major_formatter().set_useOffset(False)
        plt.grid()
        mplcursors.cursor(hover=True)
