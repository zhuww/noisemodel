import numpy as np
import scipy.linalg as sl
import os, sys
import time
from numpy import dot
from datatools.tempo import *
from fortran_utils import *
from scipy.optimize import fmin, minimize, fmin_powell

tstart = time.time()
tf = TOAfile('J1713+0747_NANOGrav_8yv0.tim')
#md = model('J1713+0747_NANOGrav_8yv0.gls.par')
#md = model('1713_noise_gls.par')
md = model('1713_noise_fmin.par')
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
p2 = [float(md.__dict__[p]) for p in T2ECORR]
np0 = len(p0)
np1 = np0 + len(p1)
np2 = np1 + len(p2)

plist = np.abs(np.array(p0 + p1 + p2))
#p0 = plist[:np0]
#p1 = plist[np0:np1]
#p2 = plist[np1:np2]

def loglikelihood(plist):
    tstart = time.time()
    """setup parameters"""
    p0 = plist[:np0]
    p1 = plist[np0:np1]
    p2 = plist[np1:np2]
    for i,p in enumerate(T2EFAC):
        md.__dict__[p] = np.abs(p0[i])
    for i,p in enumerate(T2EQUAD):
        md.__dict__[p] = np.abs(p1[i])
    for i,p in enumerate(T2ECORR):
        md.__dict__[p] = np.abs(p2[i])
    #testparfile = '1713_noise_test.par'
    #md.write(testparfile)
    #tempofit(testparfile, tf.toafile)
    #toas, freq, res, err = read_resid()
    #r = res

    """Setup S and U"""
    S_ele = []
    stack = []
    aveeph = []
    for key in md.toagrps:
        ECORR = md.__dict__['ECORR -f %s' % key]
        for epochs in md.toagrps[key]:
            S_ele.append(float(ECORR)**2)
            aveeph.append(epochs)
            idx = md.toagrps[key][epochs]
            l = np.zeros(n)
            #l[idx] = 1./len(idx)
            l[idx] = 1.
            stack.append(l)
    S_ele = np.array(S_ele)
    U = (np.vstack(stack)).T
    UT = U.T
    #print "Setup U: {0} s".format(time.time()-tstart)
    #print 'len(S_ele)', len(S_ele)

    N_ele = []
    for i, toa in enumerate(tf.toalist):
        key = ' -f %s' % (toa.flags['f']) 
        EFAC = float(md.__dict__['T2EFAC'+key])
        val = (EFAC * float(toa.TOAsigma)) ** 2
        EQUAD = float(md.__dict__['T2EQUAD'+key])
        val += EQUAD ** 2
        N_ele.append(val)
    Nvec = np.array(N_ele)

    """computing Ntilt*"""
    Nir = r / Nvec
    NiF = ((1.0/Nvec) * F.T).T
    FNiF = np.dot(F.T, NiF)
    FNir = np.dot(NiF.T, r)

    cf = sl.cho_factor(FNiF) # solve mxm system
    FNiFr = sl.cho_solve(cf, FNir)
    Ntilt_r = Nir - np.dot(NiF, FNiFr)
    d = dot(UT, Ntilt_r)

    NiU = ((1./Nvec) * U.T).T
    FNiU = np.dot(NiF.T, U)
    FNiFFNiU = sl.cho_solve(cf, FNiU)
    Ntilt_U = NiU - np.dot(NiF, FNiFFNiU)

    Nep = len(S_ele)
    Si = np.zeros((Nep, Nep))
    np.fill_diagonal(Si, 1./S_ele)
    UNtU = dot(UT , Ntilt_U)
    Sigma = (Si + UNtU)

    logdetCt = np.sum(np.log(np.array(S_ele)))

    cfSigma = sl.cho_factor(Sigma)

    logdetSigma = 2.*np.sum(np.log(np.diag(cfSigma[0])))

    logdetGNG = 2.*np.sum(np.log(np.diag(cf[0]))) + np.sum(np.log(Nvec))

    dSid = np.dot(d.T, sl.cho_solve(cfSigma, d))

    LogLike1 = 0.5 * ( dot(r.T, Ntilt_r) - dSid)
    #LogLike2 = 0.5*((n-m)*np.log(2.*np.pi) + logdetGNG + logdetSigma + logdetCt)
    LogLike2 = 0.5 * ( logdetGNG + logdetSigma + logdetCt)
    LogLike = LogLike1 + LogLike2

    #print 'r^T Ntilt^{-1} r:',dot(r.T, Ntilt_r)
    #print 'd^T Sigma^{-1] d:', dSid
    #print 'logdet(Ctilt):', logdetCt
    #print 'logdet(Sigma):', logdetSigma
    #print 'logdet(GNG):', logdetGNG
    #print 'logliklihood:', LogLike

    print "calculate likelihood: {0} s".format(time.time()-tstart) , LogLike
    print 'ECORR:', np.abs(p2)
    return LogLike

#print loglikelihood(plist)
#plist = fmin_powell(loglikelihood, plist)
plist = fmin(loglikelihood, plist)
plist = np.abs(plist)
p0 = plist[:np0]
p1 = plist[np0:np1]
p2 = plist[np1:np2]
for i,p in enumerate(T2EFAC):
    md.__dict__[p] = p0[i]
for i,p in enumerate(T2EQUAD):
    md.__dict__[p] = p1[i]
for i,p in enumerate(T2ECORR):
    md.__dict__[p] = p2[i]

md.write('1713_noise_fmin.par')
