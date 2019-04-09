import numpy as np
import coordinates, models, constants, auxiliary
from scipy import integrate
from datetime import datetime,timedelta

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as  plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import mplcursors

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
    self.energy = {"thruster": None,
                   "solar panels": None,
                   "battery": None}

class Maneuvers:

  def __init__(self,coe,spacecraft,startDate,**kwargs):

    #Perturbation Flags
    self._PERTURBATION_ATMDRAG_ = 0
    self._PERTURBATION_J2_ = 0
    self._PERTURBATION_SOLARPRESS_ = 0
    self._PERTURBATION_THRUST_ = 0
    self._PERTURBATION_MOON_ = 0
    self._PERTURBATION_SUN_ = 0
    #Other Flags
    self._INCLUDE_ENERGY_CALCULATION_ = 0

    self.spacecraft = spacecraft

    # Set up history to store data along maneuvers
    self.history = ManeuversHistory()
    self.history.r, self.history.v = coordinates.kep2cart(coe)
    self.history.coe = np.array(coe)
    self.history.mee = coordinates.kep2eq(coe)
    self.history.propMass = [self.spacecraft.wetMass-self.spacecraft.dryMass]
    self.history.datetime = [startDate]
    self.history.energy["thruster"] = [0]
    self.history.energy["solar panels"] = [0]
    self.history.energy["battery"] = [self.spacecraft.battery.energy*60*60]

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
    self.targetRun = False;
    self.targetOrbit = np.array([None, None, None, None, None])
    self.terminalConditions = {"low_altitude": True,
                               "depleted_propellant": True}

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
      print("Day:"+str(t/60/60/24)+"\tAltitude: "+str(z/1000)+" km")

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

        alpha = self.thrustProfile[0](coordinates.cart2kep(r,v))
        beta =  self.thrustProfile[1](coordinates.cart2kep(r,v))

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
  def betts(self,t,y):
    equinoctial_coordinates = y[0:6]
    propMass = y[6]
    solarPanelsEnergy = y[7]
    thrusterEnergy = y[8]
    batteryEnergy = y[9]

    p,f,g,h,k,L = equinoctial_coordinates
    
    # Inferred Parameters
    q = 1+f*np.cos(L)+g*np.sin(L)
    s = np.sqrt(1+h**2+k**2)
    r,v = coordinates.eq2cart(equinoctial_coordinates)
    z = np.linalg.norm(r) - constants.Re

    coe = coordinates.eq2kep(equinoctial_coordinates);
    #Print Day Progress
    if self.verbose:
      if t % (60*60*24) < 100:
        #print("equinocitalElements:"+str(equinoctial_coordinates))
        print("Day: {0:.3f}\tAltitude: {1:.3f} km\tMass: {2:.3f}".format(t/60/60/24, z/1000, propMass),end="\r")

    if self._INCLUDE_ENERGY_CALCULATION_:
      #Calculate Power Available
      PSolarPanels = self.spacecraft.solarPanels.power(r,self.history.datetime[0] + timedelta(seconds=t))
      if batteryEnergy > 10000:
        PBattery = self.spacecraft.battery.dischargePower
      else:
        #SIGMOID CORRECTION
        tfactor = (batteryEnergy/(10000)*2-1)*10
        PBattery = self.spacecraft.battery.dischargePower/(1+np.e**-tfactor)
        #PBattery = 0
      PAvailable = PSolarPanels + PBattery
      # Assuming other devices consume enough power to keep the battery charging very slowly
      POtherDevices = self.spacecraft.solarPanels.nominalPower*0.6
    else:
      PAvailable = 1

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
    if self._PERTURBATION_THRUST_ and propMass > 0 and PAvailable > 0: 
      if(self.targetRun):
        alpha, beta = self.getTargetAngles(coe)
      else:
        alpha = self.thrustProfile[0](coe)
        beta = self.thrustProfile[1](coe)

      mass = self.spacecraft.dryMass+propMass
      if self._INCLUDE_ENERGY_CALCULATION_:
        PThruster, thrust, massFlowRate = self.spacecraft.thruster.operationalParams(PAvailable-POtherDevices)
      else:
        PThruster = 0
        thrust = self.spacecraft.thruster.thrust
        massFlowRate = self.spacecraft.thruster.massFlowRate
      #Thrust vectoring
      RCNThrustAngle = np.array([np.cos(beta)*np.sin(alpha),
                                 np.cos(beta)*np.cos(alpha),
                                 np.sin(beta)])
      DeltaP = thrust/mass*RCNThrustAngle

      #pert = pert + Fth 
      dpropMass = -massFlowRate
    else:
      dpropMass = 0
      PThruster = 0
    
    if self._INCLUDE_ENERGY_CALCULATION_:
      # Battery Balance
      # dE = P
      dSolarPanelsEnergy = PSolarPanels
      dThrusterEnergy = PThruster
      #print(PSolarPanels,PThruster,POtherDevices)
      dBatteryEnergy = PSolarPanels-PThruster-POtherDevices
      #print("\nASDF",batteryEnergy,dBatteryEnergy)
      if (dBatteryEnergy >= 0 and batteryEnergy/60/60 >= self.spacecraft.battery.energy) or (batteryEnergy/60/60 <= 1 and dBatteryEnergy <= 0):
        #print("\nASDF",batteryEnergy,dBatteryEnergy)
        dBatteryEnergy = 0
    else:
      dSolarPanelsEnergy = 0
      dThrusterEnergy = 0
      dBatteryEnergy = 0

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
    return np.hstack((dotx,dpropMass,dSolarPanelsEnergy,dThrusterEnergy,dBatteryEnergy))

  def propagate(self,time,timestep,**kwargs):
    #Termination events functions for integrator
    def lowAltitudeBetts(t,y):
      equinoctial_coordinates = y[0:6]
      r,v = coordinates.eq2cart(equinoctial_coordinates)
      z = np.linalg.norm(r) - constants.Re
      # Zero Crossing at 200 km
      return z-200e3
    def depletedPropellantBetts(t,y):
      propMass = y[6]
      # Zero Crossing at 0
      return propMass
    #def targetReachedBetts(t,y):
    #  currentOrbit = coordinates.eq2kep(y[0:6])[0:5]
    #  elmSelector = self.targetOrbit != np.array(None)
    #  if np.any(elmSelector):
    #    currentOrbit = currentOrbit[elmSelector]
    #    targetOrbit = self.targetOrbit[elmSelector]
    #    error = abs(currentOrbit-targetOrbit)
    #    print(error)
    #    for elm in error:
    #      if elm < 1e-3:
    #        return 0
    #  return 1
      
    lowAltitudeBetts.terminal = self.terminalConditions["low_altitude"]
    depletedPropellantBetts.terminal = self.terminalConditions["depleted_propellant"]
    #targetReachedBetts.terminal = True
      
      
    print("Propagating...from day ",self.history.t[-1]/60/60/24," to ",(self.history.t[-1]+time)/60/60/24)

    t = np.linspace(self.history.t[-1],self.history.t[-1]+time,time/timestep)

    tInt = t

    if self.formulation == 'conwell':
      y0 = np.append(np.append(self.history.r[-1,:],self.history.v[-1,:]),self.history.propMass[-1])
      sol = integrate.odeint(self.conwell,y0,tInt)

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
          if perc % 10 < 1e-3:
            print("{0:.1f}%".format(perc),end="\r")
        coeLocalHist[idx,:] = coordinates.cart2kep(rLocalHist[idx],vLocalHist[idx])
      if self.verbose: print("")

    elif self.formulation == 'betts':
      y0 = np.hstack((self.history.mee[-1,:],
                      self.history.propMass[-1],
                      self.history.energy["solar panels"][-1],
                      self.history.energy["thruster"][-1],
                      self.history.energy["battery"][-1]))
      #-----INTEGRATION---------
      rtol = 1e-3
      atol = np.array([1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,
                      1e-6,
                      60*60*1e-6,60*60*1e-6,60*60*1e-6])
      events = [lowAltitudeBetts, depletedPropellantBetts]#, targetReachedBetts]
      #Using odeint
      #sol = integrate.odeint(self.betts,y0,tInt,tfirst=True,rtol=rtol,atol=atol,**kwargs)
      #Using solve_ivp
      sol = integrate.solve_ivp(self.betts,(tInt[0],tInt[-1]),y0,method="LSODA",t_eval=tInt,atol=atol,rtol=rtol,events=events,**kwargs)
      if self.verbose: print("")
      print("{0} (Status Code: {1})".format(sol.message,sol.status))
      tInt = sol.t
      sol = np.transpose(sol.y)

      y = sol[:,0:6] 
      propMass = sol[:,6]
      solarPanelsEnergy = sol[:,7]
      thrusterEnergy = sol[:,8]
      batteryEnergy = sol[:,9]


      #Initialize r and v for faster stacking
      rLocalHist = np.zeros([y.shape[0],3])
      vLocalHist = np.zeros([y.shape[0],3])
      #Initialize coe faster stacking
      coeLocalHist = np.zeros([y.shape[0],6])

      for idx,row in enumerate(y):
        if self.verbose:
          perc = idx/y.shape[0]*100
          if perc % 1 < 0.1:
            print("{0:.1f}%".format(perc),end="\r")
        r,v = coordinates.eq2cart(row)
        rLocalHist[idx,:] = r
        vLocalHist[idx,:] = v
        coeLocalHist[idx,:] = coordinates.eq2kep(row)
      if self.verbose: print("100.0%")
      print("")

      datetime = np.array([self.history.datetime[0] + timedelta(seconds=seconds) for seconds in tInt])

      # Update history
      self.history.t   = np.append(self.history.t,tInt[1:])
      self.history.coe = np.vstack((self.history.coe,coeLocalHist[1:,:]))
      self.history.r   = np.vstack((self.history.r,rLocalHist[1:,:]))
      self.history.v   = np.vstack((self.history.v,vLocalHist[1:,:]))
      self.history.datetime = np.append(self.history.datetime, datetime[1:])
      self.history.propMass = np.append(self.history.propMass,propMass[1:])
      self.history.energy["thruster"]     = np.append(self.history.energy["thruster"],thrusterEnergy[1:])
      self.history.energy["battery"]      = np.append(self.history.energy["battery"],batteryEnergy[1:])
      self.history.energy["solar panels"] = np.append(self.history.energy["solar panels"],solarPanelsEnergy[1:])

      if self.formulation == 'betts':
        self.history.mee = np.vstack((self.history.mee,y[1:,:]))

    self.history.maneuverIdxs.append(len(self.history.t))

  def impulsiveManeuver(self,dv):
    r = self.history.r[-1,:]
    v = self.history.v[-1,:]
    #dv must be in RSW Frame too
    self.history.v[-1,:] = v+coordinates.rsw2eci(r,v,dv)
    self.history.coe[-1,:] = coordinates.cart2kep(r,self.history.v[-1,:])
    self.history.mee[-1,:] = coordinates.kep2eq(self.history.coe[-1,:])

  def calculatePerturbations(self,t,r,v,propMass):
    p = [0,0,0];

    if self._PERTURBATION_ATMDRAG_:
      z = np.linalg.norm(r) - constants.Re
      #Atmospheric Drag
      vrel = v - np.cross(constants.wE,r)
      rho = models.USSA76(z)
      #rhoMin, rhoMax = models.HarrisPriester(z)
      #rho = rhoMax
      #rho = models.MSISE90(z,"mean")
      #rhoMin = models.MSISE90(z,"low")
      #rhoMax = models.MSISE90(z,"high")
      #rho = models.mixAtmosphericModels(self.history.datetime[0]+timedelta(seconds=t),rhoMin,rhoMax)
      Fd = -0.5*rho*np.linalg.norm(vrel)*vrel*(1/self.spacecraft.BC(self.spacecraft.dryMass+propMass))
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
      uhat, rS = models.solarPosition(self.history.datetime[0]+timedelta(seconds=t))
      
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

      #Radiation Pressure Coefficient (lies between 1 and 2)
      Cr = self.spacecraft.Cr
      #Spacecraft mass
      mass = self.spacecraft.dryMass+propMass
      #Radiation Pressure acceleration
      FR = -nu*PSR*Cr*As/mass*uhat

      p = p + FR

    if self._PERTURBATION_MOON_:
      r_m = models.lunarPositionAlmanac2013(self.history.datetime[0]+timedelta(seconds=t))
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
      uhat, r_s = models.solarPosition(self.history.datetime[0]+timedelta(seconds=t))
      r_sunSat = r_s-r 
      #F(q) formula from F.3 Appendix Curtis 2013
      # 	c = b - a; a << b
      # 	F = 1 - c**3/b**2
      qq = np.dot(r,(2*r_s-r))/np.linalg.norm(r_s)**2
      Fq = (qq**2 - 3*qq + 3)*qq/(1+(1-qq)**(3/2))
      pSun = constants.mu_S/np.linalg.norm(r_sunSat)**3*(Fq*r_s-r)

      p = p + pSun

    return p

  def setTargetRun(self,coe,**kwargs):
    self.targetRun = True
    self.targetOrbit = np.array(coe)
  def getTargetAngles(self,coe):
    e = coe[1]
    i = coe[2]
    omega = coe[3]
    nu = coe[5]
    E = 2*np.arctan2(np.tan(nu/2),np.sqrt((1+e)/(1-e)))
    
    alpha = np.zeros((5,))
    beta = np.zeros((5,))

    # Semi-major axis
    alpha[0] = np.arctan2(e*np.sin(nu),1+e*np.cos(nu))
    beta[0] = 0
    # Eccentricity
    alpha[1] = np.arctan2(np.sin(nu),np.cos(nu)+np.cos(E))
    beta[1] =  0
    # Inclination
    alpha[2] = 0
    beta[2] = np.sign(np.cos(omega+nu))*np.pi/2
    # Argument of Perigee
    alpha[3] = np.arctan2((1+e*np.cos(nu)),(2+e*np.cos(nu))*np.tan(nu))
    beta[3] = np.arctan2(e*np.sin(omega+nu),np.tan(i)*(np.sin(alpha[3]-nu)*(1+e*np.cos(nu))-np.cos(alpha[3])*np.sin(nu)))
    # RAAN
    alpha[4] = 0
    beta[4] = np.sign(np.sin(omega+nu))*np.pi/2
    
    Tcoe = np.zeros((3,5))
    weights= np.zeros((5,))
    for j in range(0,5):
      #Thrust vectoring
      Tcoe[:,j]= np.array([np.cos(beta[j])*np.sin(alpha[j]),
                                 np.cos(beta[j])*np.cos(alpha[j]),
                                 np.sin(beta[j])])
      if self.targetOrbit[j] == None:
        weights[j] = 0
      else:
        weights[j] = (self.targetOrbit[j]-coe[j])/abs(self.targetOrbit[j]-self.history.coe[-1,j])
        if(abs(weights[j]) < 1e-3):
          weights[j] = 0

    TOpt = np.matmul(Tcoe,weights)

    r = np.linalg.norm(TOpt)
    alphaOpt = np.arctan2(TOpt[0],TOpt[1])
    betaOpt = np.arcsin(TOpt[2]/r)

    #print("alpha,beta:\n",alpha*180/np.pi,"\n",beta*180/np.pi)
    #print("weights: ",weights)
    #print("Tcoe: \n",Tcoe)
    #print("TOpt: ",TOpt)
    #print("alphaOpt, betaOpt: ",alphaOpt*180/np.pi,betaOpt*180/np.pi)
    #print("------------------")

    return alphaOpt,betaOpt

  def calculateSecularElements(self):
    self.history.secularCoe = []
    self.history.tSecular = []
    T0 = 2*np.pi*(self.history.coe[0,0]**3/constants.mu_E)**0.5
    t = T0;
    idx0 = 0;
    idxf = (np.abs(self.history.t - t)).argmin();
    while t < self.history.t[-1]:
      #print(idx0,idxf)
      meanCoe = np.mean(self.history.coe[idx0:idxf,0:5],axis=0)
      self.history.secularCoe = np.append(self.history.secularCoe,meanCoe,axis=0);
      self.history.tSecular = np.append(self.history.tSecular,self.history.t[idx0])

      T = 2*np.pi*(self.history.coe[idxf,0]**3/constants.mu_E)**0.5
      t = t+T
      idx0 = idxf
      idxf = (np.abs(self.history.t - t)).argmin();
    self.history.secularCoe = np.reshape(self.history.secularCoe,(-1,5))

#------------------------------------------------------------------------
  def plot(self, item, itemHistory=None,**kwargs):
    useDate = False
    for key in kwargs:
      if key == "useDate":
        useDate = kwargs[key]

    if item == "coe":
      # PLOTTING CLASSICAL ORBITAL ELEMENTS
      titles = ["a","e","i","$\omega$","$\Omega$","$\\nu$"]
      ylabels = ["[km]", "", "[°]", "[°]", "[°]", "[°]"]
      timeAxis = self.history.datetime if useDate else self.history.t/60/60/24
      fig, axes = plt.subplots(3,2,figsize=(10,8),sharex=True)
      for i in range(0,6):
        for j in range(0,len(self.history.maneuverIdxs)-1):
          maneuverSlice = slice(self.history.maneuverIdxs[j],self.history.maneuverIdxs[j+1])
          if i in [2,3,4,5]:
            axes[int((i-i%2)/2),i%2].plot(timeAxis[maneuverSlice],self.history.coe[maneuverSlice,i]*180/np.pi)
          else:
            if i == 0:
              axes[int((i-i%2)/2),i%2].plot(timeAxis[maneuverSlice],self.history.coe[maneuverSlice,i]/1000)
            else:
              axes[int((i-i%2)/2),i%2].plot(timeAxis[maneuverSlice],self.history.coe[maneuverSlice,i])
          axes[int((i-i%2)/2),i%2].set_title(titles[i]+" "+ylabels[i])
          
        if useDate:
          fig.autofmt_xdate()
          axes[int((i-i%2)/2),i%2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        else:
          if i in [4,5]:
            axes[int((i-i%2)/2),i%2].set_xlabel("Tiempo [días]")
        axes[int((i-i%2)/2),i%2].yaxis.get_major_formatter().set_scientific(False)
        axes[int((i-i%2)/2),i%2].yaxis.get_major_formatter().set_useOffset(False)
        axes[int((i-i%2)/2),i%2].grid(b=True)
        if i in [0,1]:
          axes[int((i-i%2)/2),i%2].yaxis.get_major_formatter().set_scientific(True)
          axes[int((i-i%2)/2),i%2].yaxis.get_major_formatter().set_useOffset(True)
    if item == "secularCoe":
      # PLOTTING CLASSICAL ORBITAL ELEMENTS
      titles = ["a","e","i","$\omega$","$\Omega$","$\\nu$"]
      ylabels = ["[km]", "", "[°]", "[°]", "[°]", "[°]"]
      timeAxis = self.history.datetime if useDate else self.history.tSecular/60/60/24
      fig, axes = plt.subplots(3,2,figsize=(10,8),sharex=True)
      for i in range(0,5):
          if i in [2,3,4]:
            axes[int((i-i%2)/2),i%2].plot(timeAxis,self.history.secularCoe[:,i]*180/np.pi)
          else:
            if i == 0:
              axes[int((i-i%2)/2),i%2].plot(timeAxis,self.history.secularCoe[:,i]/1e3)
            else:
              axes[int((i-i%2)/2),i%2].plot(timeAxis,self.history.secularCoe[:,i])
          axes[int((i-i%2)/2),i%2].set_title(titles[i]+" "+ylabels[i])
          
          if(useDate):
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

        figXLim = ax.get_xlim()
        figYLim = ax.get_ylim()
        figZLim = ax.get_zlim()
        xx, yy = np.meshgrid(figXLim, figYLim)
        z = np.array([[0,0],[0,0]])
        ax.plot_surface(xx, yy, z, alpha=0.3,color="lightgray")

        ax.set_title("Satellite Trajectory [km]")
        ax.set_xlabel("X [km]");
        ax.set_ylabel("Y [km]");
        ax.set_zlabel("Z [km]");
      return ax

    if item == "orbitalEnergy":
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
        #fig, ax = plt.subplots(figsize=(10,4))
        #ax.plot(self.history.datetime,moonDistances);

    if item == "energyUsage":
      timeAxis = self.history.datetime if useDate else self.history.t/60/60/24
      fig, ax = plt.subplots(figsize=(10,4))
      PSolarPanels = np.diff(self.history.energy["solar panels"])/np.diff(self.history.t)
      PThruster = np.diff(self.history.energy["thruster"])/np.diff(self.history.t)
      POtherDevices = np.ones((len(self.history.t)))*(self.spacecraft.solarPanels.nominalPower)*0.6
      EBattery = self.history.energy["battery"]/60/60
      ax.plot(timeAxis[:-1], PSolarPanels,linewidth=1,label="Potencia Paneles Solares")
      ax.plot(timeAxis[:-1], PThruster,linewidth=1,label="Potencia Propulsor")
      ax.plot(timeAxis, POtherDevices,linewidth=1,label="Potencia Otros Dispositivos")
      ax.set_ylabel("Potencia [W]")
      #ax.set_ylim([-1,2])

      if useDate:
        fig.autofmt_xdate()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
      else:
        ax.set_xlabel("Tiempo [días]")
      ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
      ax.yaxis.get_major_formatter().set_scientific(False)
      ax.yaxis.get_major_formatter().set_useOffset(False)

      ax2 = ax.twinx()
      ax2.plot(timeAxis, EBattery, linewidth=1,label="Energía Baterías",color="purple")
      ax2.set_ylabel("Energía [Wh]")
      h1, l1 = ax.get_legend_handles_labels()
      h2, l2 = ax2.get_legend_handles_labels()
      ax2.legend(h1+h2,l1+l2)
      ax2.set_ylim([-2,self.spacecraft.battery.energy+2])
      ax2.grid()
      #ax2.set_ylim([-.5,.5])

    if item == "singleItem":
      timeAxis = self.history.datetime if useDate else self.history.t/60/60/24
      if np.isscalar(itemHistory) and itemHistory == None:
        raise Exception("History Data not specified.")
      else:
        fig, ax = plt.subplots(figsize=(10,4))
        for i in range(0,len(self.history.maneuverIdxs)-1):
          maneuverSlice = slice(self.history.maneuverIdxs[i],self.history.maneuverIdxs[i+1])
          ax.plot(timeAxis[maneuverSlice], itemHistory[maneuverSlice],linewidth=1)

        if useDate:
          fig.autofmt_xdate()
          ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
          ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
          ax.yaxis.get_major_formatter().set_scientific(False)
          ax.yaxis.get_major_formatter().set_useOffset(False)
        else:
          ax.set_xlabel("Time [days]")
        plt.grid()
        mplcursors.cursor(hover=True)
  def ipvPlot3D(self,**kwargs):
    import ipyvolume as ipv
    import random
    newFig = True
    for key in kwargs:
      if(key == "fig"):
        newFig = False
        fig = kwargs[key]
    if newFig:
      fig = ipv.figure(width=800, height=500);

    for i in range(0,len(self.history.maneuverIdxs)-1):
        maneuverSlice = slice(self.history.maneuverIdxs[i],self.history.maneuverIdxs[i+1])
        #randomColor = '#%02X%02X%02X' % (random.randint(0,255),random.randint(0,255),random.randint(0,255))
        colorCycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        ipv.pylab.plot(self.history.r[maneuverSlice,0]/1000, self.history.r[maneuverSlice,2]/1000, self.history.r[maneuverSlice,1]/1000,
                       color=colorCycle[i])

    if newFig:
      ipv.pylab.xlabel("X")
      ipv.pylab.ylabel("Z")
      ipv.pylab.zlabel("Y")
      ipv.pylab.style.box_off()

      #theta, phi = np.linspace(0, 2 * np.pi, 60), np.linspace(0, np.pi, 60)
      #THETA, PHI = np.meshgrid(theta, phi)
      #R = constants.Re/1e3
      #X = R * np.sin(PHI) * np.cos(THETA)
      #Y = R * np.sin(PHI) * np.sin(THETA)
      #Z = R * np.cos(PHI)

      #ipv.pylab.plot_surface(X, Z, Y, color="lightblue")
      ipv.pylab.xyzlim(10e3)
    return fig

  def ipvPlot3DMoon(self):
    import ipyvolume as ipv
    lunarPositions = np.zeros((self.history.datetime.shape[0],3))
    for idx in range(0,self.history.datetime.shape[0]):
      lunarPositions[idx,:] = models.lunarPositionAlmanac2013(self.history.datetime[idx])/1e3
    ipv.pylab.plot(lunarPositions[:,0],lunarPositions[:,2],lunarPositions[:,1])
    ipv.pylab.scatter(np.array([lunarPositions[-1,0]]),np.array([lunarPositions[-1,2]]),np.array([lunarPositions[-1,1]]),color="gray",marker="sphere")

  def calculateEclipseHours(self):
    self.history.eclipse = np.zeros((self.history.coe.shape[0],1))

    for idx in range(0,len(self.history.t)):
      r = self.history.r[idx,:]
      v = self.history.v[idx,:]
      uhat, rS = models.solarPosition(self.history.datetime[idx])
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
      self.history.eclipse[idx] = nu

  def calculatePower(self):
    self.history.solarPanelsPower = np.zeros((self.history.coe.shape[0],1))
    self.history.thrusterPower = np.zeros((self.history.coe.shape[0],1))
    for idx in range(0,len(self.history.t)):
      self.history.solarPanelsPower[idx] = self.spacecraft.solarPanels.power(self.history.r[idx,:],self.history.datetime[idx])
      self.history.thrusterPower[idx] = self.spacecraft.thruster.thrustToPower()

  def makeReport(self):
    report =  "------------MANEUVER REPORT-------------\n"
    report += "----INITIAL CONDITIONS----\n"
    report += "Date/Time:\n"
    report += "  Initial Date: "+self.history.datetime[0].strftime("%Y-%m-%d %H:%M:%S")+"\n";
    report += "\nSpacecraft:\n"
    report += "  Wet Mass: \t\t"+str(self.spacecraft.wetMass)+" kg\n"
    report += "  Dry Mass: \t\t"+str(self.spacecraft.dryMass)+" kg\n"
    report += "  Propellant Mass: \t"+str(self.spacecraft.wetMass-self.spacecraft.dryMass)+" kg\n";
    report += "  Drag Area: \t\t"+str(self.spacecraft.area)+" m2\n";
    report += "  Cd: \t\t\t"+str(self.spacecraft.Cd)+"\n";
    report += "  Cr: \t\t\t"+str(self.spacecraft.Cr)+"\n";
    report += "\nThruster:\n"
    report += "  Name/Model:\t\t"+str(self.spacecraft.thruster.name)+"\n";
    report += "  Thrust (nominal):\t"+str(self.spacecraft.thruster.thrust)+" N\n";
    report += "  Isp (nominal):\t"+str(self.spacecraft.thruster.isp)+" s\n";
    report += "  Power (nominal):\t"+str(self.spacecraft.thruster.power)+" W\n";
    report += "\nSolar Panels:\n"
    report += "  Name/Model:\t\t\t"+str(self.spacecraft.solarPanels.name)+"\n";
    report += "  Number of Panels:\t\t"+str(self.spacecraft.solarPanels.n)+"\n";
    report += "  Individual Area:\t\t"+str(self.spacecraft.solarPanels.area)+" m2\n";
    report += "  Total Power (nominal):\t"+str(self.spacecraft.solarPanels.nominalPower)+" W\n";
    report += "\nBattery:\n"
    report += "  Name/Model:\t\t"+str(self.spacecraft.battery.name)+"\n";
    report += "  Cells Configuration:\t"+str(self.spacecraft.battery.P)+"P-"+str(self.spacecraft.battery.S)+"S\n";
    report += "  Voltage:\t\t"+str(self.spacecraft.battery.voltage)+" V\n";
    report += "  Capacity:\t\t"+str(self.spacecraft.battery.capacity)+" mAh\n";
    report += "  Energy:\t\t"+str(self.spacecraft.battery.energy)+" Wh\n";
    report += "  Charge Power:\t\t"+str(self.spacecraft.battery.chargePower)+" W\n"
    report += "  Discharge Power:\t"+str(self.spacecraft.battery.dischargePower)+" W\n"
    report += "\nOrbit:\n"
    report += "  Semi-major axis (a): \t\t"+str(self.history.coe[0,0]/1e3)+" km\n"
    report += "  Eccentricity (e): \t\t"+str(round(self.history.coe[0,1],6))+"\n"
    report += "  Inclination (i): \t\t"+str(self.history.coe[0,2]*180/np.pi)+" deg\n"
    report += "  Argument of Perigee (omega): \t"+str(round(self.history.coe[0,3]*180/np.pi,5))+" deg\n"
    report += "  RAAN (Omega): \t\t"+str(round(self.history.coe[0,4]*180/np.pi,5))+" deg\n"
    report += "  True Anomaly (nu): \t\t"+str(round(self.history.coe[0,5]*180/np.pi,5))+" deg\n"
    
    
    # Report at the end of each stage
    for idx in range(1,len(self.history.maneuverIdxs)-1):
      maneuverIdx = self.history.maneuverIdxs[idx]
      report += "\n----STAGE "+str(idx)+"----\n";
      report += "Date/Time:\n"
      deltatime = self.history.datetime[maneuverIdx]-self.history.datetime[0]
      report += "  Elapsed Time:\t"+str(deltatime)+"\n";
      report += "  Date at end of stage:\t"+self.history.datetime[maneuverIdx].strftime("%Y-%m-%d %H:%M:%S")+"\n";
      report += "\nSpacecraft:\n"
      report += "  Propellant Mass: \t"+str(self.history.propMass[maneuverIdx])+" kg\n";
      report += "\nOrbit:\n"
      report += "  Semi-major axis (a): \t\t"+str(self.history.coe[maneuverIdx,0]/1e3)+" km\n"
      report += "  Eccentricity (e): \t\t"+str(round(self.history.coe[maneuverIdx,1],6))+"\n"
      report += "  Inclination (i): \t\t"+str(self.history.coe[maneuverIdx,2]*180/np.pi)+" deg\n"
      report += "  Argument of Perigee (omega): \t"+str(round(self.history.coe[maneuverIdx,3]*180/np.pi,5))+" deg\n"
      report += "  RAAN (Omega): \t\t"+str(round(self.history.coe[maneuverIdx,4]*180/np.pi,5))+" deg\n"
      report += "  True Anomaly (nu): \t\t"+str(round(self.history.coe[maneuverIdx,5]*180/np.pi,5))+" deg\n"

    report += "\n----FINAL CONDITIONS----\n";
    report += "Date/Time:\n"
    deltatime = self.history.datetime[-1]-self.history.datetime[0]
    report += "  Elapsed Time:\t"+str(deltatime)+"\n";
    report += "  End Date:\t"+self.history.datetime[-1].strftime("%Y-%m-%d %H:%M:%S")+"\n";
    report += "\nSpacecraft:\n"
    report += "  Propellant Mass: \t"+str(self.history.propMass[-1])+" kg\n";
    report += "\nOrbit:\n"
    report += "  Semi-major axis (a): \t\t"+str(self.history.coe[-1,0]/1e3)+" km\n"
    report += "  Eccentricity (e): \t\t"+str(round(self.history.coe[-1,1],6))+"\n"
    report += "  Inclination (i): \t\t"+str(self.history.coe[-1,2]*180/np.pi)+" deg\n"
    report += "  Argument of Perigee (omega): \t"+str(round(self.history.coe[-1,3]*180/np.pi,5))+" deg\n"
    report += "  RAAN (Omega): \t\t"+str(round(self.history.coe[-1,4]*180/np.pi,5))+" deg\n"
    report += "  True Anomaly (nu): \t\t"+str(round(self.history.coe[-1,5]*180/np.pi,5))+" deg\n"

    print(report)
