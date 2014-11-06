from datatools.tempo import *
from TexTable import deluxetable
from datatools.MJD import MJD_to_datetime

t = TOAfile('1713.Sep.T2.tim')
m = model('Oct.T1.RN.par')

newgroups = {}
newgroups['Mark III (L)'] = t.groups['M3A-L'] + t.groups['M3B-L']
newgroups['Mark IV (L)'] = t.groups['M4-L']
newgroups['Mark IV (S)'] = t.groups['M4-S']
newgroups['Mark IV-O (L)'] = t.groups['M4O-L']
newgroups['Mark IV-O (S)'] = t.groups['M4O-S']
newgroups['ABPP (L)'] = t.groups['ABPP-L']
newgroups['ABPP (S)'] = t.groups['ABPP-S']
newgroups['ASP (L)'] = t.groups['ASP-L']
newgroups['ASP (S)'] = t.groups['ASP-S']
newgroups['GASP (820)'] = t.groups['GASP-8']
newgroups['GASP (L)'] = t.groups['GASP-L']
newgroups['GUPPI (820)'] = t.groups['GUPPI-8']
newgroups['GUPPI (L)'] = t.groups['GUPPI-L']
newgroups['PUPPI (L)'] = t.groups['PUPPI-L']
newgroups['PUPPI (S)'] = t.groups['PUPPI-S']

oldgroups = t.groups
t.groups = newgroups
m.tempofit(t)
m.average()
 
keys = ['Mark III (L)', 'Mark IV (L)', 'Mark IV (S)', 'Mark IV-O (L)', 'Mark IV-O (S)', 'ABPP (L)','ABPP (S)','ASP (L)','ASP (S)','GASP (820)', 'GASP (L)', 'GUPPI (820)', 'GUPPI (L)', 'PUPPI (L)', 'PUPPI (S)']

BandWidth = [40, 10, 10, 10, 10, 56, 112, 64, 64, 64, 64, 800, 800, 800, 800]
Integration = [47, 58, 29, 60, 30, 60, 30, 20, 20, 20, 20, 20, 20, 20, 20]

NoOfTOA = []
for k in keys:
    NoOfTOA.append(len(t.groups[k]))

DateRange = []
def findrange(grp):
    toas = [t.toalist[i] for i in t.groups[grp]]
    mjds = [float(toa.TOA) for toa in toas]
    maxd = max(mjds)
    mind = min(mjds)
    return '%s$-$%s' % (MJD_to_datetime(mind).strftime('%Y %b'), MJD_to_datetime(maxd).strftime('%Y %b'))

for k in keys:
    krange = findrange(k)
    #print k, krange
    DateRange.append(krange)
NoOfEpoch = []
for k in keys:
    NoOfEpoch.append(len(m.averes[k]))

#print NoOfEpoch

Caption = '21 year J1713+0747 observations.'
colnames = ['System', 'Dates', 'Number of', 'Epochs', 'Bandwidth', 'Typical Integration']
data = [keys, DateRange, NoOfTOA, NoOfEpoch, BandWidth, Integration]

table = deluxetable(Caption=Caption, colnames=colnames, data=data, label='tab:obs', colsetting='lccccc')

print table
