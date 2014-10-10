import numpy as np
import scipy.linalg as sl
import os, sys
import time
from numpy import dot
from datatools.tempo import *
from fortran_utils import *
from pylab import *
from scipy.optimize import fmin, minimize, fmin_powell

secperday = 86400
dayperyear = 365.24218967
EFAC_default = 1.0
ECORR_default = 1.e-6
LARGE_NUMBER = np.exp(300)

tstart = time.time()
tf = TOAfile('1713.Sep.T2.tim')
md = model('1713_21yr_test.par')
md.tempofit(tf, DesignMatrix=True)
T2EFAC = [par for par in md.__dict__ if par.startswith('T2EFAC')]
T2EQUAD = [par for par in md.__dict__ if par.startswith('T2EQUAD')]
T2ECORR = [par for par in md.__dict__ if par.startswith('ECORR')]

if 'RNAMP' in md.__dict__:
    RNAMP = np.log(float(md.RNAMP))
else:
    RNAMP = np.log(0.06)
    md.manifest.append('RNAMP')
if 'RNIDX' in md.__dict__:
    RNIDX = float(md.RNIDX)
else:
    RNIDX = -2.17
    md.manifest.append('RNIDX')

groups = {}
for par in T2EQUAD:
    flag = ' '.join(par.split()[1:3])
    groups[flag] = tf.othergroups[flag]


md.average(lapse=0.0001, groups=groups)
t,w,A = read_design() #read in the design matrix
u, s, v = numpy.linalg.svd(A) # svd decomposition

#Np = len(s)
#F = u[..., :Np] #extract the F matrixd
#G = u[..., Np:] #extract the G matrrix

r = np.array(md.res)
n = u.shape[0]
#m = Np
M = A

p0 = [float(md.__dict__[p]) for p in T2EFAC]
p1 = [float(md.__dict__[p]) for p in T2EQUAD]
p2 = [float(md.__dict__[p]) for p in T2ECORR]
p3 = [RNAMP, RNIDX]
np0 = len(p0)
np1 = np0 + len(p1)
np2 = np1 + len(p2)
np3 = np2 + len(p3)
plist = np.array(p0 + p1 + p2 + p3)


def loglikelihood(plist):
    tstart = time.time()
    """setup parameters"""
    p0 = plist[:np0]
    p1 = plist[np0:np1]
    p2 = plist[np1:np2]
    p3 = plist[np2:np3]
    for i,p in enumerate(T2EFAC):
        md.__dict__[p] = np.abs(p0[i])
    for i,p in enumerate(T2EQUAD):
        md.__dict__[p] = np.abs(p1[i])
    for i,p in enumerate(T2ECORR):
        md.__dict__[p] = np.abs(p2[i])
    md.__dict__['RNAMP'] = p3[0]
    md.__dict__['RNIDX'] = p3[1]

    """ Setup EFAC, EQUAD"""
    N_ele = []
    EFACtags = [' '.join(tag.split()[1:3]) for tag in T2EFAC]
    EQUADtags = [' '.join(tag.split()[1:3]) for tag in T2EQUAD]
    for i, toa in enumerate(tf.toalist):
        TOAflags = [('-%s %s' % (f,toa.flags[f])) for f in toa.flags]
        try:
            key = (set(EFACtags) & set(TOAflags)).pop()
            EFAC = float(md.__dict__['T2EFAC ' + key])
        except:
            EFAC = EFAC_default
        val = (EFAC * float(toa.TOAsigma)) ** 2
        key = (set(EQUADtags) & set(TOAflags)).pop()
        EQUAD = float(md.__dict__['T2EQUAD ' + key])
        val += EQUAD ** 2
        N_ele.append(val)
    Nvec = np.array(N_ele)

    """Setup S and U for jitter parameter ECORR """
    S_ele = []
    stack = []
    aveeph = []
    for key in md.toagrps:
        try:
            ECORR = md.__dict__['ECORR %s' % key]
        except:
            ECORR = ECORR_default
        for epochs in sorted(md.toagrps[key].keys()):
            S_ele.append(float(ECORR)**2)
            aveeph.append(epochs)
            idx = md.toagrps[key][epochs]
            l = np.zeros(n)
            l[idx] = 1.
            stack.append(l)
    S_ele = np.array(S_ele)
    U = (np.vstack(stack)).T
    UT = U.T
    #print "Setup U: {0} s".format(time.time()-tstart)
    #print 'len(S_ele)', len(S_ele)

    """Setup red noise Fourier matrix, Fr"""
    N_mode = 100
    T = float(tf.end - tf.start)
    Tstart = float(tf.start)
    def RedPower(RNAMP, RNIDX, T, N_mode):
        secperyear = secperday*dayperyear
        logfac = np.log((secperyear*1e6)/(2.0*np.pi*np.sqrt(3.0)))
        T_in_s = T * secperday
        result = 2*(RNAMP-logfac) - np.log(12.) - 2.*np.log(np.pi) - np.log(T_in_s) + 3.*np.log(secperyear)
        logmu = np.log(np.array(range(1, N_mode+1))) - np.log(T_in_s ) + np.log(secperyear)
        result += logmu*RNIDX
        return np.exp(result)
    Fr_ele = np.array([[p] for p in RedPower(RNAMP, RNIDX, T, N_mode)])
    phi = np.vstack((Fr_ele, Fr_ele))
    #phi = Fr_ele

    K = np.array([ [2. * np.pi * mode / T] for mode in range(1, N_mode+1)])
    avetoa = np.hstack([md.avetoa[key]-Tstart for key in md.toagrps])
    FEsin = np.sin(avetoa * K)
    FEcos = np.cos(avetoa * K)
    FE = np.vstack((FEsin, FEcos))
    #PhiFE = phi * FE
    #C = np.dot(FE.T, PhiFE)

    """Putting the noise models together"""
    F = dot(U, FE.T) 
    #print M.shape, F.shape, U.shape
    T = np.hstack((M, F, U))
    #print T.shape

    """Putting the timing and noise parameters together """
    m = M.shape[1]
    n_f = F.shape[1]
    n_j = U.shape[1]
    n_p = T.shape[1]
    Phi_I = np.zeros((n_p, n_p))
    Pars =  [LARGE_NUMBER for i in range(m)]
    Pars += list(phi)
    Pars += list(S_ele)
    Pi = 1./np.array(Pars)
    np.fill_diagonal(Phi_I, Pi)

    """compute likelyhood"""
    d = dot(T.T, r/Nvec)
    Sigma = Phi_I + dot(T.T, (1./Nvec * T.T).T)
    #print d.shape, Sigma.shape
    cfSigma = sl.cho_factor(Sigma)

    logdetSigma = 2.*np.sum(np.log(np.diag(cfSigma[0])))
    Sid = sl.cho_solve(cfSigma, d)
    dSid = np.dot(d.T, Sid)

    logdetCt = np.sum(np.log(np.array(Pars)))
    logdetN = np.sum(np.log(Nvec))

    Nir = r / Nvec

    LogLike1 = 0.5 * ( dot(r.T, Nir) - dSid)
    #LogLike2 = 0.5*((n-m)*np.log(2.*np.pi) + logdetGNG + logdetSigma + logdetCt)
    LogLike2 = 0.5 * ( logdetN + logdetSigma + logdetCt)
    LogLike = LogLike1 + LogLike2


    print "calculate likelihood: {0} s".format(time.time()-tstart) , LogLike
    #print 'EFAC:', np.abs(p0)
    #print 'ECORR:', np.abs(p2)
    print 'RNAMP: %s, RNIDX: %s' % (np.exp(p3[0]), p3[1])

    return LogLike


#print loglikelihood(plist)
#plist = fmin(loglikelihood, plist)
plist = fmin_powell(loglikelihood, plist)
p0 = plist[:np0]
p1 = plist[np0:np1]
p2 = plist[np1:np2]
p3 = plist[np2:np3]
for i,p in enumerate(T2EFAC):
    md.__dict__[p] = np.abs(p0[i])
for i,p in enumerate(T2EQUAD):
    md.__dict__[p] = np.abs(p1[i])
for i,p in enumerate(T2ECORR):
    md.__dict__[p] = np.abs(p2[i])
md.__dict__['RNAMP'] = np.exp(p3[0])
md.__dict__['RNIDX'] = p3[1]

md.write('1713_21yr_RN.par')
