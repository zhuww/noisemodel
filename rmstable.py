from datatools.tempo import *
import numpy as np
from TexTable import *
from pylab import *
import cPickle as pickle
#from round import shortform as SF

def SF(f):
    if type(f) == str:
        return f
    else:
        return '%.3f' % f

m = model('Oct.T1.omdot.par')
t = TOAfile('1713.Sep.T2.tim')
#m.tempofit(t, GLS=True)
del m
#pickle.dump(m, open('Oct.T1.omdot.pkl', 'wb'), protocol=2)
m = pickle.load(open('Oct.T1.omdot.pkl', 'rb'))
m.groups['all'] = range(len(t.toalist))
y,err = np.load('RNrealization.npy')
#m.res -= y
#m.average()
#m.plot('date', 'averes')
#show()
#m.plot('ophase', 'averes')
#show()
print m.wrms
print m.avewrms
#pickle.dump(m, open('Oct.T1.pkl', 'wb'), protocol=2)
#keys = m.groups.keys()
#keys.sort()
sys.exit(0)
keys = ['all', 'M3A-L', 'M3B-L', 'M4-L', 'M4-S', 'M4O-L', 'M4O-S', 'ABPP-L', 'ABPP-S', 'ASP-L', 'ASP-S', 'GASP-8', 'GASP-L', 'GUPPI-8', 'GUPPI-L', 'PUPPI-L', 'PUPPI-S' ]
#print m.avewrms
wrms = [SF(m.wrms[k]) for k in keys]
avewrms = [SF(m.avewrms[k]) for k in keys]
avesig = [SF(np.mean([float(t.toalist[i].TOAsigma) for i in m.groups[k]])) for  k in keys]
EFAC = []
EQUAD = []
ECORR = []
EFACdic = {}
EQUADdic = {}
ECORRdic = {}
lines = open(m.parfile, 'r').readlines()
#print avesig
EFACline = [l.split() for l in lines if l.startswith('T2EFAC')]
EQUADline = [l.split() for l in lines if l.startswith('T2EQUAD')]
ECORRline = [l.split() for l in lines if l.startswith('ECORR')]

#print EQUADline
#sys.exit(0)

keytable = {'M3A-L':'M3A-L',
        'M3B-L':'M3B-L',
        'M4-L':'M4-L',
        'M4-S':'M4-S',
        'M4O-L':'M4O-L',
        'M4O-S':'M4O-S',
        'ABPP-L':'ABPP-L',
        'ABPP-S':'ABPP-S',
        'L-wide_PUPPI':'PUPPI-L',
        'S-wide_PUPPI':'PUPPI-S',
        'L-wide_ASP':'ASP-L',
        'S-wide_ASP':'ASP-S',
        'Rcvr1_2_GUPPI':'GUPPI-L',
        'Rcvr_800_GUPPI':'GUPPI-8',
        'Rcvr_800_GASP':'GASP-8',
        'Rcvr1_2_GASP':'GASP-L'}

for l in EFACline:
    try:
        EFACdic[keytable[l[2]]] = l[3]
    except:
        EFACdic[l[2]] = l[3]
for l in EQUADline:
    try:
        EQUADdic[keytable[l[2]]] = l[3]
    except:
        EQUADdic[l[2]] = l[3]
for l in ECORRline:
    try:
        ECORRdic[keytable[l[2]]] = l[3]
    except:
        ECORRdic[l[2]] = l[3]


for k in keys:
    if k in EFACdic:
        EFAC.append(SF(float(EFACdic[k])))
    else:
        EFAC.append('--')

    if k in EQUADdic:
        EQUAD.append(SF(float(EQUADdic[k])))
    else:
        EQUAD.append('--')

    if k in ECORRdic:
        ECORR.append(SF(float(ECORRdic[k])))
    else:
        ECORR.append('--')


WRMS = [SF(m.wrms[k]) for k in keys]
AWRMS = [SF(m.avewrms[k]) for k in keys]

savearray = (keys, WRMS, AWRMS)
np.save(open('WRMS.npy', 'wb'), savearray)

Caption = "Noise parameters and residual RMS."
colnames = ['Backends', r'$\bar{\sigma}$ ($\mu$s)', 'EFAC ($\mu$s)', 'EQUAD ($\mu$s)', 'ECORR ($\mu$s)', 'WRMS# ($\mu$s)','AWRMS# ($\mu$s)']
data = [keys, avesig, EFAC, EQUAD, ECORR, WRMS, AWRMS]
table = deluxetable(Caption=Caption, colsetting='lcccccc', colnames = colnames, data=data, label="tab:rms" , comments = ["Definition of WRMS", "(AWRMS) wegithed RMS of Epoch-averaged residual "])

print table
