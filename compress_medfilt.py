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
#threadit.func_defaults[0]['state'] = True

filename = sys.argv[1]

WINDOWSIZE = 4095

programstart = time.time()

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
tsamp = float(hdu1.header['TBIN'])
nsubint = int(hdu1.header['NAXIS2'])
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

#print obsfreq
#print type(duration)
#print type(fmin), type(fmax)

pb = PB(maxValue = 64) 

data = hdu1.data['data']
print data.shape
print data[-1,-1,-1,-1,-1]

print 'time spent reading file:', time.time() - programstart
now = time.time()

#diffdata = np.diff(data, n=WINDOWSIZE, axis=1)
#cumdata = np.cumsum(data, axis=1) 
#print cumdata.shape
#print 'time spent sum the data:', time.time() - now
#diffdata = np.diff(cumdata, n=WINDOWSIZE, axis=1)
#print 'time spent diff the data:', time.time() - now

l,m,n,o,p = data.shape
data = data.reshape(l,m,n,o/4,4,p).sum(axis=4).astype(np.float32)
nchan = o/4
print 'time spent quenching the frequency:', time.time() - now
now = time.time()

polX = []
polY = []

def func(i):
    #dataX = data[i,:,0,:,0].reshape((-1,nchan))
    #dataY = data[i,:,1,:,0].reshape((-1,nchan))
    dataX = data[:,:,0,i,0].flatten()
    dataY = data[:,:,1,i,0].flatten()
    #print dataX.shape, dataY.shape
    #datapolX = medfilt2d(dataX, (WINDOWSIZE,1))
    #datapolY = medfilt2d(dataY, (WINDOWSIZE,1))
    datapolX = medfilt(dataX, WINDOWSIZE)
    datapolY = medfilt(dataY, WINDOWSIZE)
    #polX.append(datapolX)
    #polY.append(datapolY)
    pb(i)
    return datapolX + datapolY

#for i in range(64):

    #dataX = data[i,:,0,:,0].reshape((-1,4096))
    #dataY = data[i,:,1,:,0].reshape((-1,4096))
    #print dataX.shape, dataY.shape
    #datapolX = medfilt2d(dataX, (WINDOWSIZE,1))
    #datapolY = medfilt2d(dataY, (WINDOWSIZE,1))
    #polX.append(datapolX)
    #polY.append(datapolY)
    #pb(i)


#res = np.array(threadit(func, [[i] for i in range(nchan)])).T
res = np.array(threadit(func, [[i] for i in range(nchan)])).T
print 'time spent medfilt the data:', time.time() - now

print res.shape
#res = res.reshape((2,-1))

#datapolX = res[:,0,:,:].squeeze()
#datapolY = res[:,1,:,:].squeeze()
dataX = data[:,:,0,:,:].squeeze()
dataY = data[:,:,1,:,:].squeeze()
dataS = (dataX + dataY).reshape((-1, nchan))

onebitdata = (res - dataS).astype(np.uint8)

print onebitdata.shape

#onebitX = (dataX > datapolX).astype(np.uint8)
#onebitY = (dataY > datapolY).astype(np.uint8)
        

#data = data.sum(axis=2)
#data = data.squeeze()
#l,m,n = data.shape
#data = data.reshape(l*m,n)
#m,n = data.shape
#data = data.reshape(m/64, 64, n)
#data = data.mean(axis=1)
#print data.shape
#print data[0,0,0,0]

#print onebitX.shape
#data = onebitX.reshape(-1,4096,nchan)
#l,m,n = data.shape
#data = data.reshape(l*m,n)
#m,n = data.shape
#data = data.reshape(m/64, 64, n)
#data = data.mean(axis=1)

#'''

data = data.reshape(-1,4096,nchan)
l,m,n = data.shape
data = data.reshape(l*m,n)
m,n = data.shape
data = data.reshape(m/64, 64, n)
data = data.mean(axis=1)

from pylab import *
imshow(data.T, aspect='auto', origin='low')  #, extent=(0., float(duration), float(fmin), float(fmax)))
xlabel('t (s)')
ylabel('frequency (MHz)')
show()


data = onebitdata.reshape(-1,4096,nchan)
l,m,n = data.shape
data = data.reshape(l*m,n)
m,n = data.shape
data = data.reshape(m/64, 64, n)
data = data.mean(axis=1)

#'''
#from pylab import *
imshow(data.T, aspect='auto', origin='low')  #, extent=(0., float(duration), float(fmin), float(fmax)))
#xlabel('t (s)')
#ylabel('frequency (MHz)')
show()
#'''
