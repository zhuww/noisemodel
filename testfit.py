from pylab import *
from datatools.tempo import *

m = model('1713.Sep.T2.par')
t = TOAfile('1713.Sep.T2.tim')
#m.thawall()
m.tempo2fit(t)
m.average()
print m.chisq, m.dof

a1 = subplot(211)
a2 = subplot(212)


m.plot('date', 'DMX', ax=a1)
m.plot('date', 'averes', ax=a2)
#m.plot('mjd', 'DMX')
#m1 = model(m.newpar.parfile)
#del m,t
#m1.plot('mjd', 'DMX')
show()
