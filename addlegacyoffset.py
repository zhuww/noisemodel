from datatools.tempo import *

tf = TOAfile('1713.Jan.T2.tim')

for toa in tf.toalist:
    #if toa.flags.has_key('label') and toa.flags['label'] == 'Legacy':
        #try:
            #print toa.flags['to'] 
        #except:
            #print toa.flags['i']
    if toa.flags.has_key('i') and toa.flags['i'] in ['ABPP-L', 'M4-L', 'M4O-L']:
        if toa.flags.has_key('to'):
            toa.flags['to'] = str(Decimal(toa.flags['to']) + Decimal('-0.000000049'))
        else: 
            toa.flags['to'] = '-0.000000049'

print tf.tempo2fmt()
