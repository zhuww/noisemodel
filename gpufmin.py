import numpy as np
import scipy.linalg as sl
from cholesky import Gcholesky, Gchosolve
import os, sys
import time
from numpy import dot
from datatools.tempo import *
from fortran_utils import *
from pylab import *
from scipy.optimize import fmin, fmin_powell
from rankreduced import get_rr_rep, pl_psd
np.set_printoptions(precision=3, suppress=True)

secperday = 86400
dayperyear = 365.24218967
secperyear = secperday*dayperyear
EFAC_default = 1.0
ECORR_default = 1.e-6
EQUAD_default = 1.e-6
LARGE_NUMBER = np.exp(50)
SMALL_NUMBER = np.exp(-100)
ampfac = ((secperyear*1e6)/(2.0*np.pi*np.sqrt(3.0)))

tstart = time.time()
tf = TOAfile('1713.Sep.T2.tim')
Tstart = float(tf.start)
Tspan = float(tf.end - tf.start)/dayperyear #in unit of year
toas = np.array([(float(toa.TOA)-Tstart)/dayperyear for toa in tf.toalist]) #in unit of year



md = model('1713_21yr_test.par')
#md = model('1713.Oct.test.par')
#md = model('1713_21yr_omdot.par')
#md = model('1713_21yr_JAE.par')
md.tempofit(tf, DesignMatrix=True)
T2EFAC = [par for par in md.__dict__ if par.startswith('T2EFAC')]
T2EQUAD = [par for par in md.__dict__ if par.startswith('T2EQUAD')]
T2ECORR = [par for par in md.__dict__ if par.startswith('ECORR')]

if 'RNAMP' in md.__dict__:
    RNAMP = np.abs(float(md.RNAMP))
else:
    RNAMP = np.abs(0.06)
    md.manifest.append('RNAMP')
if 'RNIDX' in md.__dict__:
    RNIDX = float(md.RNIDX)
else:
    RNIDX = -2.17
    md.manifest.append('RNIDX')

groups = {}
flaglist = []
for par in T2EQUAD:
    flag = ' '.join(par.split()[1:3])
    groups[flag] = tf.othergroups[flag]
    flaglist.append(flag)


md.average(lapse=0.0001, groups=groups)
t,w,M = read_design() #read in the design matrix
#u, s, v = numpy.linalg.svd(A) # svd decomposition

#Np = len(s)
#F = u[..., :Np] #extract the F matrixd
#G = u[..., Np:] #extract the G matrrix

r = np.array(md.res)
n = r.shape[0]
#m = Np
#M = A


""" Setup matrix U"""
ECORRtags = [' '.join(tag.split()[1:3]) for tag in T2ECORR]
S_idx = []
stack = []
aveeph = []
ECORRlist = []
for key in T2ECORR:
    flag = ' '.join(key.split()[1:3])
    if key  in md.__dict__:
        ECORRlist.append(float(md.__dict__[key]))
        S_idx.append(0)
        for epochs in sorted(md.toagrps[flag].keys()):
            aveeph.append(epochs)
            idx = md.toagrps[flag][epochs]
            S_idx[-1] += 1
            l = np.zeros(n)
            l[idx] = 1.
            stack.append(l)
        S_idx[-1] = np.array(S_idx[-1])
S_ele = np.hstack([ ECORRlist[i]**2 *np.ones(k) for i,k in enumerate(S_idx)])
U = (np.vstack(stack)).T
UT = U.T

""" Setup EFAC, EQUAD"""
EFACtags = [' '.join(tag.split()[1:3]) for tag in T2EFAC]
EQUADtags = [' '.join(tag.split()[1:3]) for tag in T2EQUAD]
EFACidx = {}
EQUADidx = {}
SIGMA = []
EFAClist = []
EQUADlist = []
for tag in EFACtags:
    EFACidx[tag] = []
    EFAC = float(md.__dict__['T2EFAC ' + tag])
    EFAClist.append(EFAC)
for tag in EQUADtags:
    EQUADidx[tag] = []
    EQUAD = float(md.__dict__['T2EQUAD ' + tag])
    EQUADlist.append(EQUAD)
for i, toa in enumerate(tf.toalist):
    SIGMA.append(float(toa.TOAsigma))
    TOAflags = [('-%s %s' % (f,toa.flags[f])) for f in toa.flags]
    try:
        key = (set(EFACtags) & set(TOAflags)).pop()
        EFACidx[key].append(i)
    except:
        EFAC = EFAC_default
    key = (set(EQUADtags) & set(TOAflags)).pop()
    EQUADidx[key].append(i)

"""Setup red noise Fourier matrix, Fr"""
f, F = get_rr_rep(toas, Tspan, 1./Tspan/4.7, 50, 20)
#fyr = 1./secperyear


p0 = [float(md.__dict__[p]) for p in T2EFAC]
p1 = [np.log(float(md.__dict__[p])) for p in T2EQUAD]
p2 = [np.log(float(md.__dict__[p])) for p in T2ECORR]
p3 = [np.log(RNAMP), RNIDX]
np0 = len(p0)
np1 = np0 + len(p1)
np2 = np1 + len(p2)
np3 = np2 + len(p3)
plist = np.array(p0 + p1 + p2 + p3)


def loglikelihood(plist, logspace=True):
    tstart = time.time()
    """setup parameters"""
    p0 = np.array(plist[:np0])
    p1 = np.array(plist[np0:np1])
    p2 = np.array(plist[np1:np2])
    p3 = np.array(plist[np2:np3])
    if logspace:
        CoordTransTerm =  (np.sum(p1) + np.sum(p2) + p3[0]) #sume all the logterm for coordinate transformation dz = z d (ln(z))
        p1 = np.exp(p1)
        p2 = np.exp(p2)
        p3[0] = np.exp(p3[0])
    else:
        CoordTransTerm = 0.
    """Load parameters"""
    efac = np.ones(n) 
    for i,tag in enumerate(EFACtags):
        efac[EFACidx[tag]] = p0[i]
    Nvec = np.array(SIGMA * efac)**2
    equad = np.zeros(n)
    for i,tag in enumerate(EQUADtags):
        equad[EQUADidx[tag]] = p1[i]
    Nvec += (equad**2)
    S_ele = np.hstack([ p2[i]**2. *np.ones(k) for i,k in enumerate(S_idx)])
    RNAMP = np.abs(p3[0]) #/ ampfac
    RNIDX = p3[1]
    #phi = (RNAMP**2)/12/np.pi/np.pi *fyr**(-3-RNIDX) * f**RNIDX
    phi = (RNAMP**2) * f**RNIDX #assuming f is is 1/yr units.

    """Putting the noise models together"""
    T = np.hstack((M, F, U))

    """Putting the timing and noise parameters together """
    m = M.shape[1]
    n_f = F.shape[1]
    n_j = U.shape[1]
    n_p = T.shape[1]
    Phi_I = np.zeros((n_p, n_p))
    Pars =  [LARGE_NUMBER for i in range(m)]
    Pars += list(phi)
    Pars += list(S_ele)
    Pi = 1./np.array(Pars) #+ SMALL_NUMBER
    np.fill_diagonal(Phi_I, Pi)

    """compute likelyhood"""
    d = dot(T.T, r/Nvec)
    Sigma = Phi_I + dot(T.T, (1./Nvec * T.T).T)
    #print d.shape, Sigma.shape
    try:
        #cfSigma = sl.cholesky(Sigma)
        #np.save('sigma', Sigma)
        cfSigma = Gcholesky(Sigma)
        #print "pass cholesky"
    #except LinAlgError:
    except :
        print Sigma



    logdetSigma = 2.*np.sum(np.log(np.diag(cfSigma)))
    #Sid = sl.cho_solve((cfSigma, False), d)
    Sid = Gchosolve(cfSigma, d)
    #print "pass solve"
    dSid = np.dot(d.T, Sid)

    logdetCt = np.sum(np.log(np.array(Pars[m:])))
    logdetN = np.sum(np.log(Nvec))

    Nir = r / Nvec

    LogLike1 = 0.5 * ( dot(r.T, Nir) - dSid)
    #LogLike2 = 0.5*((n-m)*np.log(2.*np.pi) + logdetGNG + logdetSigma + logdetCt)
    LogLike2 = 0.5 * ( logdetN + logdetSigma + logdetCt)
    #Penalty= Gamma * np.max(plist**4)
    #CoordTransTerm =  0 #(np.sum(p1) + np.sum(p2)) #sume all the logterm for coordinate transformation dz = z d (ln(z))
    #LogLike = LogLike1 + LogLike2 + Penalty
    LogLike = LogLike1 + LogLike2 - CoordTransTerm


    print "calculate likelihood: {0} s".format(time.time()-tstart) , LogLike, p3
    #print 'EFAC:', np.abs(p0)
    #print 'ECORR:', np.abs(p2)
    #print 'EQUAD:', np.abs(p1)
    #print 'RNAMP: %s, RNIDX: %s' % (np.exp(p3[0]), p3[1])
    #print 'plist:', LogLike, CoordTransTerm, (p3)
    #print 'phi', np.log(phi)
    #print 'f', np.log10(f)

    return LogLike

#LogLike = lambda x:loglikelihood(x) * -1.
#from testmcmc import SliceSampleMC as mcmc
#res, xmax, pmax = mcmc(LogLike, plist, m=100, n = 10000, ss=0.5,progressbar=True)
#res = np.array(res)
#np.save('rednoise', res)
#plist = xmax
#sys.exit(0)
plist = fmin(loglikelihood, plist)
#plist = fmin_powell(loglikelihood, plist)
p0 = plist[:np0]
p1 = plist[np0:np1]
p2 = plist[np1:np2]
p3 = plist[np2:np3]
for i,p in enumerate(T2EFAC):
    md.__dict__[p] = np.abs(p0[i])
for i,p in enumerate(T2EQUAD):
    md.__dict__[p] = np.exp(p1[i])
for i,p in enumerate(T2ECORR):
    md.__dict__[p] = np.exp(p2[i])
md.__dict__['RNAMP'] = np.exp(p3[0])
md.__dict__['RNIDX'] = p3[1]

md.write('1713_21yr_gpu.par')
