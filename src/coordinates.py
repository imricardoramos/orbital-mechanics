import numpy as np
import constants

def kep2cart(coe):
  """ Keplerian to Cartesian """
  a = coe[0]
  e = coe[1]
  i = coe[2]
  omega = coe[3]
  Omega = coe[4]
  nu = coe[5]

  E = 2*np.arctan2(np.tan(nu/2),np.sqrt((1+e)/(1-e)))

  mu = constants.mu_E
  rc = a*(1-e*np.cos(E))
  o = rc*np.array([np.cos(nu),np.sin(nu)])
  odot = (mu*a)**0.5/rc*np.array([-np.sin(E), (1-e**2)**0.5*np.cos(E)])
  rotationMatrix =  np.array(
          [[np.cos(omega)*np.cos(Omega)-np.sin(omega)*np.cos(i)*np.sin(Omega), -(np.sin(omega)*np.cos(Omega)+np.cos(omega)*np.cos(i)*np.sin(Omega))],
           [np.cos(omega)*np.sin(Omega)+np.sin(omega)*np.cos(i)*np.cos(Omega),   np.cos(omega)*np.cos(i)*np.cos(Omega)-np.sin(omega)*np.sin(Omega)],
           [np.sin(omega)*np.sin(i), np.cos(omega)*np.sin(i)]])
  r = np.matmul(rotationMatrix,o)
  v = np.matmul(rotationMatrix,odot)

  return r,v

def cart2kep(r,v):
  """ Cartesian to Keplerian """
  mu = constants.mu_E
  h = np.cross(r,v)
  rNorm = np.linalg.norm(r)
  vNorm = np.linalg.norm(v)
  eVect = np.cross(v,h)/mu-r/rNorm
  e = np.linalg.norm(eVect)
  n = np.array([-h[1],h[0],0]);
  nu = np.arccos(np.dot(eVect,r)/(e*rNorm))
  if np.dot(r,v) < 0:
    nu = 2*np.pi-nu
  i = np.arccos(h[2]/np.linalg.norm(h))
  E = 2*np.arctan2(np.tan(nu/2),np.sqrt((1+e)/(1-e)))
  Omega = np.arccos(n[0]/np.linalg.norm(n))
  if n[1] < 0:
    Omega = 2*np.pi-Omega
  omega = np.arccos(np.dot(n,eVect)/(np.linalg.norm(n)*e))
  if eVect[2] < 0:
    omega = 2*np.pi-omega
  #M = E-e*np.sin(E)
  a = 1/(2/rNorm-vNorm**2/mu)

  # Ranges:
  # i = [0,180°]
  # omega (ARGP) = [0,360°]
  # Omega (RAAN) = [0,360°]
  # M = [0,360°]

  return np.array([a,e,i,omega,Omega,nu])

def kep2eq(coe):
  """ Keplerian to Equinoctial """
  a = coe[0]
  e = coe[1]
  i = coe[2]
  omega = coe[3]
  Omega = coe[4]
  nu = coe[5]

  p = a*(1-e**2)
  f = e*np.cos(omega+Omega)
  g = e*np.sin(omega+Omega)
  h = np.tan(i/2)*np.cos(Omega)
  k = np.tan(i/2)*np.sin(Omega)
  L = omega+Omega+nu

  return np.array([p,f,g,h,k,L])

def eq2cart(mee):
  """ Equinoctial to Cartesian """
  p = mee[0]
  f = mee[1]
  g = mee[2]
  h = mee[3]
  k = mee[4]
  L = mee[5]

  alpha2 = h**2-k**2
  s2 = 1+h**2+k**2
  w = 1+f*np.cos(L)+g*np.sin(L)
  r = p/w

  r = np.array([(r/s2)*(np.cos(L)+alpha2*np.cos(L)+2*h*k*np.sin(L)),
                (r/s2)*(np.sin(L)-alpha2*np.sin(L)+2*h*k*np.cos(L)),
                (2*r/s2)*(h*np.sin(L)-k*np.cos(L))])
  v = np.array([-1/s2*np.sqrt(constants.mu_E/p)*((1+alpha2)*(np.sin(L)+g)-2*h*k*(np.cos(L)+f)),
                -1/s2*np.sqrt(constants.mu_E/p)*((alpha2-1)*(np.cos(L)+f)+2*h*k*(np.sin(L)+g)),
                2/s2*np.sqrt(constants.mu_E/p)*(h*(np.cos(L)+f)+k*(np.sin(L)+g))])

  return r,v

def eq2kep(mee):
  p = mee[0]
  f = mee[1]
  g = mee[2]
  h = mee[3]
  k = mee[4]
  L = mee[5]

  a = p/(1-f**2-g**2)
  e = np.sqrt(f**2+g**2)
  i = np.arctan2(2*np.sqrt(k**2+h**2),1-h**2-k**2)
  omega = np.arctan2(g*h-f*k,f*h+g*k)
  Omega = np.arctan2(k,h)
  nu = L-Omega-omega

  omega = (omega + 2*np.pi)/np.pi % 2*np.pi
  Omega = (Omega + 2*np.pi)/np.pi % 2*np.pi
  nu = nu/np.pi % 2*np.pi

  return np.array([a,e,i,omega,Omega,nu])

def rsw2eci(r,v,vector):
  #RSW frame definition from r and v
  rhat = r/np.linalg.norm(r)
  w = np.cross(rhat,v)
  what = w/np.linalg.norm(w)
  shat = np.cross(what,rhat)
  rsw = [rhat, shat, what]
  return np.dot(np.transpose(rsw),vector)
def eci2rsw(r,v,vector):
  #RSW frame definition from r and v
  rhat = r/np.linalg.norm(r)
  w = np.cross(rhat,v)
  what = w/np.linalg.norm(w)
  shat = np.cross(what,rhat)
  rsw = [rhat, shat, what]
  return np.dot(rsw,vector)
