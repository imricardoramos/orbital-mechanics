import numpy as np
from coordinates import kep2cart, kep2eq, eq2cart, eq2kep, cart2kep
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
    self.history.propMass = [self.spacecraft.wetMass-self.spacecraft.dryMass]
    self.history.datetime = startDate

    self.history.r.shape = (1,3)
    self.history.v.shape = (1,3)
    self.history.coe.shape = (1,6)
    self.history.mee.shape = (1,6)

    # Other options
    self.verbose = False;
    self.formulation = "betts";
    for key in kwargs:
      if key == 'verbose':
        self.verbose = kwargs[key]
      if key == 'formulation':
        self.formulation = kwargs[key]

  def addPerturbation(self,perturbation):
    if perturbation == "atmosphere":
      self._PERTURBATION_ATMDRAG_ = 1
    if perturbation == "J2":
      self._PERTURBATION_J2_ = 1
    if perturbation == "solar_pressure":
      self._PERTURBATION_SOLARPRESS_ = 1
    if perturbation == "thrust":
      self._PERTURBATION_THRUST_ = 1
      def a(coe):
        e = coe[1]
        nu = coe[5]
        alpha = np.arctan2(e*np.sin(nu),(1+e*np.cos(nu)))
        return alpha;
      def b(coe):
        return 0;
      self.thrustProfile = (a,b)
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
    propMass = stateVector[6]
    z = np.linalg.norm(r)-constants.Re

    if t % (60*60*24) < 100 and self.verbose:
      print("Day:"+str(t/60/60/24)+"\tHeight: "+str(z/1000)+" km")

    if(z > 100e3):
      #Calculate Perturbations
      p = self.calculatePerturbations(t,r,v,propMass)


      if self._PERTURBATION_THRUST_ and propMass > 0:
        #ECI frame definition from r and v
        rhat = r/np.linalg.norm(r)
        w = np.cross(rhat,v)
        what = w/np.linalg.norm(w)

        that = v/np.linalg.norm(v)
        r2 = np.cross(what,that)
        r2hat = r2/np.linalg.norm(r2)

        rtn = np.array([r2hat, that, what])

        alpha = self.thrustProfile[0](cart2kep(r,v))
        beta =  self.thrustProfile[1](cart2kep(r,v))

        mass = self.spacecraft.dryMass+propMass
        
        #Thrust vectoring
        RCNThrustAngle = np.array([np.cos(beta)*np.sin(alpha),
                                   np.cos(beta)*np.cos(alpha),
                                   np.sin(beta)])
        Fth = np.dot(self.spacecraft.thruster.thrust/mass*RCNThrustAngle,rtn)
        p = p + Fth 
        dpropMass = -self.spacecraft.thruster.massFlowRate
      else:
        dpropMass = 0
      
      E = np.linalg.norm(v)**2/2-constants.mu_E/np.linalg.norm(r)
      if t % (60*60*24) < 100 and self.verbose:
        print("Orbital Energy:"+str(E));
      #Differential Equations
      drdt = v
      dvdt = -constants.mu_E*(r/(np.linalg.norm(r)**3))+p

      return np.append(np.append(drdt,dvdt),dpropMass)
    else:
      return np.append(np.append(-r,-v),dpropMass)
  def betts(self,y,t):
    equinoctial_coordinates = y[0:6]
    p = equinoctial_coordinates[0]
    f = equinoctial_coordinates[1]
    g = equinoctial_coordinates[2]
    h = equinoctial_coordinates[3]
    k = equinoctial_coordinates[4]
    L = equinoctial_coordinates[5]

    propMass = y[6]
    
    # Infer Parameters
    q = 1+f*np.cos(L)+g*np.sin(L)
    s = np.sqrt(1+h**2+k**2)
    r,v = eq2cart(equinoctial_coordinates)
    z = np.linalg.norm(r) - constants.Re
    
    #Print Day Progress
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

    #Calculate Perturbations
    pert = self.calculatePerturbations(t,r,v,propMass)
    DeltaP = 0;

    # THRUST PERTURBATON HAD TO BE CALCULATED HERE...
    if self._PERTURBATION_THRUST_ and propMass > 0:
      #ECI frame definition from r and v
      #that = v/np.linalg.norm(v)
      #r2 = np.cross(what,that)
      #r2hat = r2/np.linalg.norm(r2)

      #rtn = [r2hat, that, nhat]

      alpha = self.thrustProfile[0](eq2kep(equinoctial_coordinates))
      beta = self.thrustProfile[1](eq2kep(equinoctial_coordinates))

      mass = self.spacecraft.dryMass+propMass
      
      #Thrust vectoring
      RCNThrustAngle = np.array([np.cos(beta)*np.sin(alpha),
                                 np.cos(beta)*np.cos(alpha),
                                 np.sin(beta)])
      DeltaP = self.spacecraft.thruster.thrust/mass*RCNThrustAngle

      #pert = pert + Fth 
      dpropMass = -self.spacecraft.thruster.massFlowRate
    else:
      dpropMass = 0

    #Transform perturbations to RSW Frame
    Delta = np.dot(rsw,pert) + DeltaP

    A = np.array([[0, 2*p/q*np.sqrt(p/constants.mu_E),0],
                  [np.sqrt(p/constants.mu_E)*np.sin(L), np.sqrt(p/constants.mu_E)*1/q*((q+1)*np.cos(L)+f), -np.sqrt(p/constants.mu_E)*g/q*(h*np.sin(L)-k*np.cos(L))],
                  [-np.sqrt(p/constants.mu_E)*np.cos(L), np.sqrt(p/constants.mu_E)*(1/q)*((q+1)*np.sin(L)+g), np.sqrt(p/constants.mu_E)*(f/q)*(h*np.sin(L)-k*np.cos(L))],
                  [0, 0, np.sqrt(p/constants.mu_E)*s**2*np.cos(L)/(2*q)],
                  [0, 0, np.sqrt(p/constants.mu_E)*s**2*np.sin(L)/(2*q)],
                  [0, 0, np.sqrt(p/constants.mu_E)*(1/q)*(h*np.sin(L)-k*np.cos(L))]])
    b = np.array([0, 0, 0, 0, 0, np.sqrt(constants.mu_E*p)*(q/p)**2])

    dotx = np.matmul(A,Delta) + b
    return np.append(dotx,dpropMass)

  def propagate(self,time,timestep):
    print("Propagating...from day ",self.history.t[-1]/60/60/24," to ",(self.history.t[-1]+time)/60/60/24)

    #Integrate
    t = np.linspace(self.history.t[-1],self.history.t[-1]+time,time/timestep)

    if self.formulation == 'conwell':
      y0 = np.append(np.append(self.history.r[-1,:],self.history.v[-1,:]),self.history.propMass[-1])
      sol = integrate.odeint(self.conwell,y0,t)

      #Initialize r and v for faster stacking
      rLocalHist = np.zeros([sol.shape[0],3])
      vLocalHist = np.zeros([sol.shape[0],3])
      #Initialize coe faster stacking
      coeLocalHist = np.zeros([sol.shape[0],6])
      
      rLocalHist = sol[:,0:3]
      vLocalHist = sol[:,3:6]
      propMass = sol[:,6]
      for idx,row in enumerate(sol[:,0:6]):
        if self.verbose:
          perc = idx/sol.shape[0]*100
          if perc % 10 == 0:
            print(str(perc)+"%")
        coeLocalHist[idx,:] = cart2kep(rLocalHist[idx],vLocalHist[idx])

    elif self.formulation == 'betts':
      y0 = np.append(self.history.mee[-1,:],self.history.propMass[-1])
      sol = integrate.odeint(self.betts,y0,t)
      y = sol[:,0:6] 
      propMass = sol[:,6]

      #Initialize r and v for faster stacking
      rLocalHist = np.zeros([y.shape[0],3])
      vLocalHist = np.zeros([y.shape[0],3])
      #Initialize coe faster stacking
      coeLocalHist = np.zeros([y.shape[0],6])

      for idx,row in enumerate(y):
        if self.verbose:
          perc = idx/y.shape[0]*100
          if perc % 10 == 0:
            print(str(perc)+"%")
        r,v = eq2cart(row)
        rLocalHist[idx,:] = r
        vLocalHist[idx,:] = v
        coeLocalHist[idx,:] = eq2kep(row)

    datetime = np.array([self.startDate + timedelta(seconds=seconds) for seconds in t])

    # Update history
    self.history.t   = np.append(self.history.t,t)
    if self.formulation == 'betts':
      self.history.mee = np.vstack((self.history.mee,y))
    self.history.coe = np.vstack((self.history.coe,coeLocalHist))
    self.history.r   = np.vstack((self.history.r,rLocalHist))
    self.history.v   = np.vstack((self.history.v,vLocalHist))
    self.history.maneuverIdxs.append(len(self.history.t))
    self.history.propMass = np.append(self.history.propMass,propMass)
    self.history.datetime = np.append(self.history.datetime,datetime)

#  def propagate3(self,time):
#    print("Propagating...from day ",self.history.t[-1]/60/60/24," to ",(self.history.t[-1]+time)/60/60/24)
#    #Integrate
#    y0 = np.append(self.history.mee[-1,:],self.history.propMass[-1])
#    t = np.linspace(self.history.t[-1],self.history.t[-1]+time,time/60)
#    sol = integrate.solve_ivp(self.betts,(t[0],t[-1]),y0,method="RK45")
#    t = sol.t
#    y = sol.y[:,0:6] 
#    propMass = sol.y[:,6]
#
#    #Initialize r and v for faster stacking
#    rLocalHist = np.zeros([y.shape[0],3])
#    vLocalHist = np.zeros([y.shape[0],3])
#
#    for idx,row in enumerate(y):
#      perc = idx/y.shape[0]*100
#      if perc % 10 == 0:
#        print(str(perc)+"%")
#      r,v = eq2cart(row)
#      rLocalHist[idx,:] = r
#      vLocalHist[idx,:] = v
#    datetime = np.array([self.startDate + timedelta(seconds=seconds) for seconds in t])
#
#    # Update history
#    self.history.t   = np.append(self.history.t,t)
#    self.history.mee = np.vstack((self.history.mee,y))
#    self.history.r   = np.vstack((self.history.r,rLocalHist))
#    self.history.v   = np.vstack((self.history.v,vLocalHist))
#    self.history.maneuverIdxs.append(len(self.history.t))
#    self.history.propMass = np.append(self.history.propMass,propMass)
#    self.history.datetime = np.append(self.history.datetime,datetime)

  def impulsiveManeuver(self,dv):
    r = self.history.r[-1,:]
    v = self.history.v[-1,:]
    #RSW frame definition from r and v
    rhat = r/np.linalg.norm(r)
    w = np.cross(rhat,v)
    what = w/np.linalg.norm(w)
    shat = np.cross(what,rhat)
    rsw = [rhat, shat, what]
    #dv must be in RSW Frame too
    self.history.v[-1,:] = v+np.dot(rsw,dv)
    self.history.coe[-1,:] = cart2kep(r,self.history.v[-1,:])
    self.history.mee[-1,:] = kep2eq(self.history.coe[-1,:])

  def calculatePerturbations(self,t,r,v,propMass):
    p = [0,0,0];

    if self._PERTURBATION_ATMDRAG_:
      z = np.linalg.norm(r) - constants.Re
      #Atmospheric Drag
      vrel = v - np.cross(constants.wE,r)
      Fd = -0.5*models.atmosDensity(z/1000)*np.linalg.norm(vrel)*vrel*(1/self.spacecraft.BC(self.spacecraft.dryMass+propMass))
      p = p + Fd

    if self._PERTURBATION_J2_:
      J2 = 0.00108263

      zz = r[2]

      rNorm = np.linalg.norm(r)
      bigTerm = np.array([1/rNorm*(5*zz**2/rNorm**2-1),
                          1/rNorm*(5*zz**2/rNorm**2-1),
                          1/rNorm*(5*zz**2/rNorm**2-3)])

      #J2 Gravity Gradient
      FJ2 = (3/2)*J2*constants.mu_E*constants.Re**2/rNorm**4*bigTerm*r
      p = p + FJ2

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
      CR = 2
      #Spacecraft mass
      mass = self.spacecraft.dryMass+propMass
      #Radiation Pressure acceleration
      FR = -nu*PSR*CR*As/mass*uhat

      p = p + FR

    if self._PERTURBATION_MOON_:
      r_m = models.lunarPositionAlmanac2013(self.startDate+timedelta(seconds=t))
      r_ms = r_m-r 
      #pMoon = constants.mu_M*(r_ms/np.linalg.norm(r_ms)**3-r_m/np.linalg.norm(r_m)**3)

      #F(q) formula from F.3 Appendix Curtis 2013
      # 	c = b - a; a << b
      # 	F = 1 - c**3/b**2
      qq = np.dot(r,(2*r_m-r))/np.linalg.norm(r_m)**2
      Fq = (qq**2 - 3*qq + 3)*qq/(1+(1-qq)**(3/2))
      pMoon = constants.mu_M/np.linalg.norm(r_ms)**3*(Fq*r_m-r)

      p = p + pMoon 

    if self._PERTURBATION_SUN_:
      uhat, r_s = models.solarPosition(self.startDate+timedelta(seconds=t))
      r_sunSat = r_s-r 
      #F(q) formula from F.3 Appendix Curtis 2013
      # 	c = b - a; a << b
      # 	F = 1 - c**3/b**2
      qq = np.dot(r,(2*r_s-r))/np.linalg.norm(r_s)**2
      Fq = (qq**2 - 3*qq + 3)*qq/(1+(1-qq)**(3/2))
      pSun = constants.mu_S/np.linalg.norm(r_sunSat)**3*(Fq*r_s-r)

      p = p + pSun


    return p
#------------------------------------------------------------------------
  def plot(self, item, itemHistory=None,**kwargs):
    if item == "coe":
      # PLOTTING CLASSICAL ORBITAL ELEMENTS
      titles = ["a","e","i","$\omega$","$\Omega$","$\\nu$"]
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
      for key in kwargs:
        if key == "ax":
          axExtern = kwargs["ax"]
      #Plot 3D Trajectory
      if 'axExtern' not in locals():
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
      else:
        ax = axExtern

      markers = np.zeros([len(self.history.maneuverIdxs)-1,3])
      
      for i in range(0,len(self.history.maneuverIdxs)-1):
          maneuverSlice = slice(self.history.maneuverIdxs[i],self.history.maneuverIdxs[i+1])
          ax.plot3D(self.history.r[maneuverSlice,0]/1000,
                    self.history.r[maneuverSlice,1]/1000,
                    self.history.r[maneuverSlice,2]/1000,linewidth=1)
          markers[i,:]= self.history.r[self.history.maneuverIdxs[i],:]/1000
      ax.plot3D(markers[:,0],markers[:,1],markers[:,2],"k.")

      if 'axExtern' not in locals():
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
      return ax

    if item == "energy":
        fig, ax = plt.subplots(figsize=(10,4))
        moonDistances = np.array([]);
        for num,date in enumerate(self.history.datetime):
          moonVector = models.lunarPositionAlmanac2013(date)
          moonDistances = np.append(moonDistances,np.linalg.norm(moonVector-self.history.r[num]))
          
        earthEnergy = np.linalg.norm(self.history.v,axis=1)**2/2 - constants.mu_E/np.linalg.norm(self.history.r,axis=1)
        moonEnergy = np.linalg.norm(self.history.v,axis=1)**2/2 - constants.mu_M/moonDistances
        ax.plot(self.history.datetime,earthEnergy,label="Earth Energy")
        ax.plot(self.history.datetime,moonEnergy, label="Moon Energy")
        fig.legend()
        fig.autofmt_xdate()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        plt.grid()
        fig, ax = plt.subplots(figsize=(10,4))
        ax.plot(self.history.datetime,moonDistances);

    if item == "singleItem":
      if np.isscalar(itemHistory) and itemHistory == None:
        raise Exception("History Data not specified.")
      else:
        fig, ax = plt.subplots(figsize=(10,4))
        for i in range(0,len(self.history.maneuverIdxs)-1):
            maneuverSlice = slice(self.history.maneuverIdxs[i],self.history.maneuverIdxs[i+1])
            ax.plot(self.history.datetime[maneuverSlice], itemHistory[maneuverSlice],linewidth=1)

        fig.autofmt_xdate()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.yaxis.get_major_formatter().set_scientific(False)
        ax.yaxis.get_major_formatter().set_useOffset(False)
        plt.grid()
        mplcursors.cursor(hover=True)
  def setTargetOrbit(self,coe_f):
    T = 0;
    for idx in len(coe):
      coe_0 = self.coe[0:idx]
      coe_cur = self.coe[-1:idx]
      T = T + (1-deltaK)*(coe_f[idx]-coe_cur)/(coe_0-coe_cur)*T_coe
