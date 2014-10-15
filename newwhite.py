import numpy as np
import scipy.linalg as sl
import os, sys
import time
from numpy import dot
from datatools.tempo import *
from fortran_utils import *
from pylab import *
from scipy.optimize import fmin, minimize, fmin_powell
from rankreduced import get_rr_rep, pl_psd

secperday = 86400
dayperyear = 365.24218967
secperyear = secperday*dayperyear
EFAC_default = 1.0
ECORR_default = 1.e-6
LARGE_NUMBER = np.exp(300)

tstart = time.time()
tf = TOAfile('1713.Sep.T2.tim')
Tstart = float(tf.start)
Tspan = float(tf.end - tf.start)
toas = np.array([float(toa.TOA)-Tstart for toa in tf.toalist])



#md = model('1713_21yr_test.par')
#md = model('1713_21yr_JAE.par')
md = model('1713.Sep.NM.par')
md.tempofit(tf, DesignMatrix=True)
T2EFAC = [par for par in md.__dict__ if par.startswith('T2EFAC')]
T2EQUAD = [par for par in md.__dict__ if par.startswith('T2EQUAD')]
T2ECORR = [par for par in md.__dict__ if par.startswith('ECORR')]

#if 'RNAMP' in md.__dict__:
    #RNAMP = np.log(float(md.RNAMP))
#else:
    #RNAMP = np.log(0.06)
    #md.manifest.append('RNAMP')
#if 'RNIDX' in md.__dict__:
    #RNIDX = float(md.RNIDX)
#else:
    #RNIDX = -2.17
    #md.manifest.append('RNIDX')

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


""" Setup matrix U"""
S_ele = []
stack = []
aveeph = []
for key in md.toagrps:
    try:
        ECORR = md.__dict__['ECORR %s' % key]
        for epochs in sorted(md.toagrps[key].keys()):
            S_ele.append(float(ECORR)**2)
            aveeph.append(epochs)
            idx = md.toagrps[key][epochs]
            l = np.zeros(n)
            l[idx] = 1.
            stack.append(l)
    except:
        ECORR = ECORR_default
        #do nothing
S_ele = np.array(S_ele)
U = (np.vstack(stack)).T
UT = U.T




p0 = [float(md.__dict__[p]) for p in T2EFAC]
p1 = [float(md.__dict__[p]) for p in T2EQUAD]
p2 = [float(md.__dict__[p]) for p in T2ECORR]
#p3 = [RNAMP, RNIDX]
np0 = len(p0)
np1 = np0 + len(p1)
np2 = np1 + len(p2)
#np3 = np2 + len(p3)
#plist = np.array(p0 + p1 + p2 + p3)
plist = np.array(p0 + p1 + p2 )


def loglikelihood(plist, Gamma=1.e-5):
    tstart = time.time()
    """setup parameters"""
    p0 = plist[:np0]
    p1 = plist[np0:np1]
    p2 = plist[np1:np2]
    #p3 = plist[np2:np3]
    for i,p in enumerate(T2EFAC):
        md.__dict__[p] = np.abs(p0[i])
    for i,p in enumerate(T2EQUAD):
        md.__dict__[p] = np.abs(p1[i])
    for i,p in enumerate(T2ECORR):
        md.__dict__[p] = np.abs(p2[i])
    #md.__dict__['RNAMP'] = np.exp(p3[0])
    #md.__dict__['RNIDX'] = p3[1]

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
    for key in md.toagrps:
        try:
            ECORR = md.__dict__['ECORR %s' % key]
            S_ele += [float(ECORR)**2] * len(md.toagrps[key])
        except:pass
        #try:
            #ECORR = md.__dict__['ECORR %s' % key]
            #for epochs in sorted(md.toagrps[key].keys()):
                #S_ele.append(float(ECORR)**2)
        #except:
            #ECORR = ECORR_default
    S_ele = np.array(S_ele)
    #print "Setup U: {0} s".format(time.time()-tstart)
    #print 'len(S_ele)', len(S_ele)

    #f, F = get_rr_rep(toas, Tspan, 1./Tspan/4.7, 50, 20)
    #fyr = 1./secperyear
    #phi = np.exp(RNAMP*2)/12/np.pi/np.pi *fyr**(RNIDX+3) * f**RNIDX

    """Putting the noise models together"""
    #T = np.hstack((M, F, U))
    T = np.hstack((M, U))

    """Putting the timing and noise parameters together """
    m = M.shape[1]
    #n_f = F.shape[1]
    n_j = U.shape[1]
    n_p = T.shape[1]
    Pars =  [LARGE_NUMBER] * m 
    #Pars += list(phi)
    Pars += list(S_ele)
    Phi_I = np.zeros((n_p, n_p))
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

    logdetCt = np.sum(np.log(np.array(Pars[m:])))
    logdetN = np.sum(np.log(Nvec))

    Nir = r / Nvec

    LogLike1 = 0.5 * ( dot(r.T, Nir) - dSid)
    #LogLike2 = 0.5*((n-m)*np.log(2.*np.pi) + logdetGNG + logdetSigma + logdetCt)
    LogLike2 = 0.5 * ( logdetN + logdetSigma + logdetCt)
    Penalty= Gamma * np.max(plist**4)
    LogLike = LogLike1 + LogLike2 + Penalty


    print "calculate likelihood: {0} s".format(time.time()-tstart) , LogLike
    #print 'EFAC:', np.abs(p0)
    print 'ECORR:', np.abs(p2)
    #print 'RNAMP: %s, RNIDX: %s' % (np.exp(p3[0]), p3[1])

    return LogLike


LogLike = lambda x:loglikelihood(x) * -1.
from testmcmc import SliceSampleMC as mcmc
from pylab import *
res = mcmc(LogLike, plist)
scatter(res[:,-1], res[:,-2])
show()

sys.exit(0)
#print loglikelihood(plist)
plist = fmin(loglikelihood, plist)
#plist = fmin_powell(loglikelihood, plist)

p0 = plist[:np0]
p1 = plist[np0:np1]
p2 = plist[np1:np2]
#p3 = plist[np2:np3]
for i,p in enumerate(T2EFAC):
    md.__dict__[p] = np.abs(p0[i])
for i,p in enumerate(T2EQUAD):
    md.__dict__[p] = np.abs(p1[i])
for i,p in enumerate(T2ECORR):
    md.__dict__[p] = np.abs(p2[i])
#md.__dict__['RNAMP'] = np.exp(p3[0])
#md.__dict__['RNIDX'] = p3[1]

md.write('1713_21yr_WN.par')
