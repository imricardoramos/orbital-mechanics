import numpy as np
from constants import constants

def kep2cart(com):
    """ Keplerian to Cartesian """
    a = com[0]
    e = com[1]
    i = com[2]
    omega = com[3]
    Omega = com[4]
    M0 = com[5]
    t = com[6]
    mu = constants.mu_E

    t0=0
    deltat = 86400*(t-t0)
    M = M0+deltat*(mu/a**3)**0.5
    
    for E in range(0,3600):
        E = E/10*np.pi/180
        if(abs(E-e*np.sin(E)-M) < 0.01):
            break;
    
    nu = 2*np.arctan2((1+e)**0.5*np.sin(E/2),(1-e)**0.5*np.cos(E/2))
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
    a = 1/(2/rNorm-vNorm**2/mu)

    return a,e

def kep2eq(a,e,i,omega,Omega,nu):
  """ Keplerian to Equinoctial """
  p = a*(1-e**2)
  f = e*np.cos(omega+Omega)
  g = e*np.sin(omega+Omega)
  h = np.tan(i/2)*np.cos(omega)
  k = np.tan(i/2)*np.sin(omega)
  L = omega+Omega+nu

  return [p,f,g,h,k,L]

def eq2cart(p,f,g,h,k,L):
  """ Equinoctial to Cartesian """
  alpha = np.sqrt(h**2-k**2)
  s = np.sqrt(1+h**2+k**2)
  w = 1+f*np.cos(L)+g*np.sin(L)
  r = p/w
  x = (r/s**2)*(np.cos(L)+alpha**2*np.cos(L)+2*h*k*np.sin(L))
  y = (r/s**2)*(np.sin(L)+alpha**2*np.sin(L)+2*h*k*np.cos(L))
  z = (r/s**2)*(h*np.sin(L)-k*np.cos(L))

  return [x,y,z]
