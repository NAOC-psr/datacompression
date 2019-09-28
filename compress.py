import numpy as np 
#import pyfits
import astropy.io.fits as pyfits
import os
#import datetime
import time
import sys
from scipy.signal import medfilt2d, medfilt
from ProgressBar import progressBar as PB
from threadit import spamit as threadit
from collections import deque
#threadit.func_defaults[0]['state'] = True

filename = sys.argv[1]

WINDOWSIZE = 4096 

print """
Running mean onebitise into filterbank file with a window size of %d samples.
""" % WINDOWSIZE

now = time.time()

hdulist = pyfits.open(filename)
hdu0 = hdulist[0]
hdu1 = hdulist[1]
#data1 = hdu1.data['data']
header1 = hdu1.header

obsfreq = hdu0.header['OBSFREQ']
obsnchan = hdu0.header['OBSNCHAN']
obsbw = hdu0.header['OBSBW']
fmin = float(obsfreq - obsbw/2.)
fmax = float(obsfreq + obsbw/2.)
nf = obsnchan
df = hdu1.header['CHAN_BW']
tsamp = float(hdu1.header['TBIN'])
nsubint = int(hdu1.header['NAXIS2'])
samppersubint = int(hdu1.header['NSBLK'])
nsamp = nsubint * samppersubint
sourcename = hdu0.header['SRC_NAME']
ra = hdu0.header['RA']
dec = hdu0.header['DEC']
RA = ''.join(ra.split(':')).replace(' ', '')
DEC = ''.join(dec.split(':')).replace(' ', '')
src_raj = float(RA)
src_dej = float(DEC)
tstart = float("%d.%d" % (hdu0.header['STT_IMJD'], hdu0.header['STT_SMJD']))
nbits = hdu0.header['BITPIX']
obsmjd = float(str(int(hdu0.header['STT_IMJD'])) + '.' + str(int(hdu0.header['STT_SMJD'])))
Azimuth = hdu1.data['TEL_AZ'][0]
Zenith = hdu1.data['TEL_ZEN'][0]
header = hdu0.header + hdu1.header
dtype = ''
duration = tsamp * nsamp
fch1 = hdu1.data['DAT_FREQ'][0][-1]
data = hdu1.data['data'].sum(axis=2)

data = data.reshape((-1, nf))
nsamp, nchan = data.shape


from runningmean import runningmean
import ctypes

print 'time took to setting up.', time.time() - now

dataout = np.zeros(data.shape, dtype='u1')

now = time.time()
dataout = runningmean(data, WINDOWSIZE)
print 'time took to running mean:', time.time() -now

#l,m = dataout.shape
#dataout = dataout.reshape(l/64,64,m).sum(axis=1)
#from pylab import *
#imshow(dataout.T, aspect='auto')
#show()

now = time.time()
data = np.packbits(dataout.reshape(-1,8)[:,::-1])
print 'bitpacking took:', time.time() - now

now = time.time()

from filwriter import filterbank_saver
outname = filename[:-5] + '.fil'

fout = open(outname, 'wb')
ibeam = 0
nbeam = 1
nbits = 1
df *= -1
filterbank_saver(fout, sourcename, ibeam, nbeam, obsmjd, tstart, Azimuth, Zenith, src_raj, src_dej, nchan, nbits, tsamp, fch1, df)
fout.write(data)
fout.close()
print 'save file took:', time.time() - now


