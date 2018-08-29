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
  n = float(line2Elms[7])*2*np.pi/(60*60*24)
  a = (constants.mu_E/n**2)**(1/3)
  coe = [a,e,i,omega,Omega,M]
  if int(line1Elms[3][0:2]) < 60:
    year = int("20"+line1Elms[3][0:2])
  else:
    year = int("19"+line1Elms[3][0:2]) 
  date = datetime.datetime(year,1,1)+datetime.timedelta(days=float(line1Elms[3][2:]))

  return coe,date
