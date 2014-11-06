import numpy as np
from pylab import *
import pickle

md = pickle.load(open('Oct.T1.pkl', 'rb'))
md.average()
#print md.wrms
#print md.avewrms

#sys.exit(0)

idx = md.groups['PUPPI-S']
res = md.res[idx]
err = md.err[idx]
wgt = 1./err**2
wg2 = md.weight[idx]
dwg = wgt - wg2
wgt = wg2
#print np.mean(wgt), np.mean(wg2), np.mean(dwg)
#print np.allclose(wgt, wg2)
#print wg2
wmres = np.sum(res*wgt)/np.sum(wgt)
mres = np.mean(res)
sumres = np.sum(res)
wsum = np.sum(wgt)
#wrms = sqrt((sum(res**2*wgt) - wmres*sumres)/wsum)
wrms = sqrt(sum(res**2*wgt)/wsum - wmres**2)
print wrms
