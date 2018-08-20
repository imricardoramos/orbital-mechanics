import numpy as np
from constants import constants
def kep2cart(com):
    a = com[0]
    e = com[1]
    i = com[2]
    omega = com[3]
    Omega = com[4]
    M0 = com[5]
    t = com[6]
    mu = com[7]

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
    h = np.cross(r,v)
    rNorm = np.linalg.norm(r)
    vNorm = np.linalg.norm(v)
    eVect = np.cross(v,h)/constants.mu_E-r/rNorm
    e = np.linalg.norm(eVect)
    a = 1/(2/rNorm-vNorm**2/constants.mu_E)
    return a,e
