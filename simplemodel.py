import numpy as np
import scipy.linalg as sl
import os, sys
import time
from numpy import dot
from datatools.tempo import *
from fortran_utils import *
from scipy.optimize import fmin, fmin_powell

tstart = time.time()
#tf = TOAfile('J1713+0747_NANOGrav_8yv0.tim')
#md = model('J1713+0747_NANOGrav_8yv0.gls.par')
tf = TOAfile('1713.Jul.T2.tim')
md = model('1713.Jul.T1.par')
#md = model('1713_efac_equad.par')
md.tempofit(tf, DesignMatrix=True)
#md.average(lapse=0.0001)
md.average()
t,w,A = read_design() #read in the design matrix
u, s, v = numpy.linalg.svd(A) # svd decomposition
#u,s,v = np.load('DesignMatrix_u.npy')
Np = len(s)
F = u[..., :Np] #extract the F matrixd
G = u[..., Np:] #extract the G matrrix

r = np.array(md.res)
n = u.shape[0]
m = Np
M = A
#Nvec = 1./w

T2EFAC = [par for par in md.__dict__ if par.startswith('T2EFAC')]
T2EQUAD = [par for par in md.__dict__ if par.startswith('T2EQUAD')]
T2ECORR = [par for par in md.__dict__ if par.startswith('ECORR')]

p0 = [float(md.__dict__[p]) for p in T2EFAC]
p1 = [float(md.__dict__[p]) for p in T2EQUAD]
#p2 = [float(md.__dict__[p]) for p in T2ECORR]
np0 = len(p0)
np1 = np0 + len(p1)
#np2 = np1 + len(p2)

plist = p0 + p1 #+ p2

def loglikelihood(plist):
    tstart = time.time()
    """setup parameters"""
    p0 = plist[:np0]
    p1 = plist[np0:np1]
    #p2 = plist[np1:np2]
    for i,p in enumerate(T2EFAC):
        md.__dict__[p] = p0[i]
    for i,p in enumerate(T2EQUAD):
        md.__dict__[p] = p1[i]
    #for i,p in enumerate(T2ECORR):
        #md.__dict__[p] = p2[i]
    

    N_ele = []
    for i, toa in enumerate(tf.toalist):
        key = ' -f %s' % (toa.flags['f']) 
        EFAC = float(md.__dict__['T2EFAC'+key])
        val = (EFAC * float(toa.TOAsigma)) ** 2
        EQUAD = float(md.__dict__['T2EQUAD'+key])
        val += EQUAD ** 2
        N_ele.append(val)
    Nvec = np.array(N_ele)

    Nir = r / Nvec
    rNir = dot(r.T, Nir)
    chisq = 0.5 * rNir
    logdetN = 0.5 * np.sum(np.log(Nvec))
    print 'chisq:', chisq, 'logdetN', logdetN
    NegLogLike = chisq + logdetN

    print 'EFAC:', p0
    return NegLogLike

#print loglikelihood(plist)
plist = fmin_powell(loglikelihood, plist)
p0 = plist[:np0]
p1 = plist[np0:np1]
#p2 = [0] * 8
for i,p in enumerate(T2EFAC):
    md.__dict__[p] = p0[i]
for i,p in enumerate(T2EQUAD):
    md.__dict__[p] = p1[i]
#for i,p in enumerate(T2ECORR):
    #md.__dict__[p] = p2[i]
md.write('1713_efac_equad.par')
