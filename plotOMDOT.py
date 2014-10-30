from pylab import *
from datatools.tempo import *

m = model('1713.Oct.omdot.par')
t = TOAfile('1713.Sep.T2.tim')
m.tempofit(t, GLS=True)
m.average()
m.plot('mjd', 'averes')
show()
