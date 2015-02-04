from pylab import *
from datatools.tempo import *

#m = model('1713.Sep.T2.par')
#m = model('1713_21yr_JSE.par')
#m = model('Oct.T1.OMDOT.par')
#m = model('Jan.T1.par')
m = model('1713.Feb.T2.par')
t = TOAfile('1713.Feb.T2.tim')
#m = model('overlap_L.par')
#t = TOAfile('overlap_L.tim')
#m.thawall()
m.tempo2fit(t) #, GLS=True)
m.average()
print m.chisq, m.dof

a1 = subplot(211)
a2 = subplot(212)

#m.plot('date', 'DMX', ax=a1)
m.plot('mjd', 'averes', ax=a2, LegendOn=False)
#m.plot('mjd', 'DMX')
m1 = model(m.newpar.parfile)
del m,t
m = model('1713.Oct.T2.par')
t = TOAfile('1713.Sep.T2.tim')
m.tempo2fit(t) #, GLS=True)
m.average()
print m.chisq, m.dof
m.plot('mjd', 'averes', ax=a1, LegendOn=False)
#m1.plot('mjd', 'DMX', ax=a1)
show()
