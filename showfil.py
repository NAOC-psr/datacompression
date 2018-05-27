from sigpyproc.Readers import FilReader as filterbank
import numpy as np
import sys, os
from pylab import *

fil = filterbank(sys.argv[1])

fch1 = fil.header['fch1']
df = fil.header['foff']
fmin = fil.header['fbottom']
fmax = fil.header['ftop']
nsamp = fil.header['nsamples']
tsamp = fil.header['tsamp']
nf = fil.header['nchans']
tstart = fil.header['tstart']
nchans = fil.header['nchans']
hdrlen = fil.header['hdrlen']

print 'fch1:', fch1
print 'df:', df
print 'fmin:', fmin
print 'fmax:', fmax
print 'nsamp:', nsamp
print 'tsamp:', tsamp
print 'nf:', nf
print 'nchans:', nchans
print 'tstart:', tstart

print 'number of channels:', nf
fil._file.seek(hdrlen)
data = fil._file.cread(nchans*nsamp)
data = np.array(data.reshape((nsamp, nf)).transpose(), order='C')


l, m = data.shape
data = data.reshape( (l, m/64, 64) ).sum(axis=2)
print data.shape

imshow(data, aspect='auto')
show()
