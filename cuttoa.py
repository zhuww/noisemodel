from datatools.tempo import *

tf = TOAfile('1713.Jan.T2.tim')

tf.start = Decimal(53342)
tf.end = Decimal(53454)

with open('overlap_L.tim', 'wt') as f:
    f.write(tf.tempo2fmt(applycut=True))

del tf
tf = TOAfile('all_reprocessed_ASP.toa')
with open('ASP_1713_ol.toa', 'wt') as f:
    for toa in tf.toalist:
        if toa.flags['tmplt'] == 'J1713+0747.L-wide.PUPPI.9y.x.sum.sm' and toa.TOA < 53393:
            f.write(toa.tempo2fmt() + '\n')
