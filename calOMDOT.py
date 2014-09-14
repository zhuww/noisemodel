from datatools.tempo import *
import numpy as np
secperday = 24*3600 
dayperyear = 365.24218967
pi = np.pi
#G = 6.673e-11
#Msun = 1.98892e30
#c = 2.99792458e8
twopi = 6.283185307179586
#fac = 1.536e-16 
Tsun = 4.925490947e-6


#pf = model('1713.Dec.mcmc.par')
#pf = model('1713.May.OMDOT.par')
#pf = model('1713.May.par')
#pf = model('1713.Jun.OMDOT.par')
pf = model('J1713+0747.par')
#pf = model('1713.9yr.S.par')
#pf = model('J1713+0747_NANOGrav_8yv0.gls.par')
#pf = model('1713.GLS.out.par')

M2 = float(pf.M2[0])
SINI = float(pf.SINI[0])
Pb = float(pf.PB[0]) * secperday
a = float(pf.A1[0]) 
#M1 = (Pb/2/pi*sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M1 = Pb/2/pi*((Tsun*(M2*SINI)**3/a**3)**(0.5))-M2
#M1 = 1.4
E = float(pf.E[0])
f = 4*pi**2*a**3/Pb**2/Tsun
print 'Timing M1,M2,Mtot,f:',M1,M2, M1+M2, f #, (M2*SINI)**3/(M1+M2)**2
print 'SINI:', SINI
print 'Timing OMDOT:', pf.OMDOT[0], '(',pf.OMDOT[1],')', 'degress per year'

#print M1+M2, E, SINI

OMDOT = 3.*(2*pi/Pb)**(5./3)*(Tsun*(M1+M2))**(2./3)/(1. - E**2)
OMDOT_PA = OMDOT * 180./pi * secperday * dayperyear
print 'periastron advance:', OMDOT_PA  , 'degress per year'

Omega = float(str(pf.PAASCNODE))
thetamu = np.arctan(float(str(pf.PMRA[0]/pf.PMDEC[0]))) #in radian
#thetamu = 180. + np.arctan(float(str(pf.PMRA[0]/pf.PMDEC[0])))/np.pi*180
#print 'thetamu:', thetamu, 'Omega:',Omega
mu = np.sqrt(float(pf.PMRA[0])**2 + float(pf.PMDEC[0])**2)
muprime = mu/float(pf.SINI[0])
OMDOT_SM = muprime*np.cos((thetamu-Omega*pi/180.))/60/60 /1000 

print 'Thetamu - Omega:', (thetamu-Omega*pi/180.)*180./pi

print 'change of viewing angle:', OMDOT_SM,'degrees per year'

print 'Observed OMDOT:',  pf.OMDOT[0], 'Expected OMDOT:', OMDOT_PA+OMDOT_SM

#print 3.*(2*pi/(secperday*88))**(5./3)*(1.)*(2./3)/(1. - E**2)* 180./pi * secperday * 365 * 60

print 'Test code using Mercury'
M1 = 1.
M2 = 3.3e23/1.989e30
Pb = 87.969 * secperday
E = 0.205630

OMDOT = 3.*(2*pi/Pb)**(5./3)*(Tsun*(M1+M2))**(2./3)/(1. - E**2)
print 'OMDOT per century:', OMDOT * 100* 180./pi * secperday * dayperyear * 60 * 60, 'seconds'