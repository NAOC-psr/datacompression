import numpy as np 
#import pyfits
import astropy.io.fits as pyfits
import os
import datetime
import time
import sys

filename = sys.argv[1]

hdulist = pyfits.open(filename)
#print len(hdulist) 
hdu0 = hdulist[0]
hdu1 = hdulist[1]
#data1 = hdu1.data['data']
header1 = hdu1.header
#print hdu0.header
#print hdu1.header


#print hdu0.header['OBSFREQ'], hdu0.header['OBSNCHAN'], hdu0.header['OBSBW']
#print hdu1.header['NCHAN'], hdu1.header['CHAN_BW'], hdu1.header['TBIN'], hdu1.header['NBITS']
#print data1.shape
#print hdulist['SUBINT'].data[0]['DAT_FREQ']
#print hdulist['SUBINT'][:12]

#fchannel = hdulist['SUBINT'].data[0]['DAT_FREQ']
#fch1 = fchannel[0]
obsfreq = hdu0.header['OBSFREQ']
obsnchan = hdu0.header['OBSNCHAN']
obsbw = hdu0.header['OBSBW']
fmin = float(obsfreq - obsbw/2.)
fmax = float(obsfreq + obsbw/2.)
nf = obsnchan
df = hdu1.header['CHAN_BW']
tsamp = hdu1.header['TBIN']
nsubint = hdu1.header['NAXIS2']
samppersubint = int(hdu1.header['NSBLK'])
nsamp = nsubint * samppersubint
sourename = hdu0.header['SRC_NAME']
ra = hdu0.header['RA']
dec = hdu0.header['DEC']
tstart = float("%d.%d" % (hdu0.header['STT_IMJD'], hdu0.header['STT_SMJD']))
nbits = hdu0.header['BITPIX']
header = hdu0.header + hdu1.header
dtype = ''
duration = tsamp * nsamp

print obsfreq
print type(duration)
print type(fmin), type(fmax)

data = hdu1.data['data']
data = data.sum(axis=2)
data = data.squeeze()
l,m,n = data.shape
data = data.reshape(l*m,n)
m,n = data.shape
data = data.reshape(m/64, 64, n)
data = data.mean(axis=1)
print data.shape

from pylab import imshow, show, xlabel, ylabel, plot
#imshow(data.T, aspect='auto', origin='low', extent=(0., float(duration), float(fmin), float(fmax)))
plot(data.sum(axis=1))
xlabel('t (s)')
ylabel('frequency (MHz)')
show()
