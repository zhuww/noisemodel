from datatools.tempo import *
import numpy as np
import scipy.optimize as opt
from pylab import *

pf = PARfile('J1713+0747.par')

twopi = 2*np.pi
fac = 1.536e-16 * 1.e12
x = float(str(pf.A1[0]))
sini = float(str(pf.SINI[0]))
cosi = -1. * np.sqrt(1 - sini**2)
mu = np.sqrt(float(str(pf.PMRA[0]**2+pf.PMDEC[0]**2)))
thetamu = 180. + np.arctan(float(str(pf.PMRA[0]/pf.PMDEC[0])))/np.pi*180
xdot = float(str(pf.XDOT[0]))
fxdot = lambda Omega: -1.* fac * x * mu * (cosi/sini) * sin((thetamu-Omega)*twopi/360.) - xdot
print 'mu:', mu, 'thetamu:', thetamu, 
print 'xdot:', xdot, 'fxdot:', fxdot(92.5)
omega = np.arange(0, 360, 30)
plot(omega, fxdot(omega))
show()
print opt.brentq(fxdot, 70, 130)
