import numpy as np
from datatools.tempo import *
from pylab import *

md = model('overlap_L.par')
tf = TOAfile('overlap_L.tim')

md.freezeall('SINI')
md.freezeall('M2')
md.freezeall('XDOT')
md.tempofit(tf)
md.average()

a1 = subplot(211)
a2 = subplot(212)

md.plot('mjd', 'DMX', ax=a1)
md.plot('mjd', 'averes', LegendOn=True, ax=a2, LegendLoc=4)
show()
