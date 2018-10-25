from constants import constants
import numpy as np
import datetime
def parseTle(tlefile):
  with open(tlefile,'r') as f:
      name = f.readline()[:-1]
      line1 = f.readline()[:-1]
      line2 = f.readline()
  line1Elms = list(filter(None,line1.split(" ")))
  line2Elms = list(filter(None,line2.split(" ")))

  i = float(line2Elms[2])*np.pi/180
  Omega = float(line2Elms[3])*np.pi/180
  e = float("0."+line2Elms[4])
  omega = float(line2Elms[5])*np.pi/180
  M = float(line2Elms[6])*np.pi/180

  for E in range(0,3600):
      E = E/10*np.pi/180
      if(abs(E-e*np.sin(E)-M) < 0.01):
          break;
  
  nu = 2*np.arctan2((1+e)**0.5*np.sin(E/2),(1-e)**0.5*np.cos(E/2))

  n = float(line2Elms[7])*2*np.pi/(60*60*24)
  a = (constants.mu_E/n**2)**(1/3)

  if int(line1Elms[3][0:2]) < 60:
    year = int("20"+line1Elms[3][0:2])
  else:
    year = int("19"+line1Elms[3][0:2]) 

  coe = [a,e,i,omega,Omega,nu]
  date = datetime.datetime(year,1,1)+datetime.timedelta(days=float(line1Elms[3][2:]))

  return coe,date

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
