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

programstart = time.time()

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
RA = ''.join(ra.split(':'))
DEC = ''.join(dec.split(':'))
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

pb = PB(maxValue = 64) 

data = hdu1.data['data']

print 'time spent reading file:', time.time() - programstart
now = time.time()

data = data.sum(axis=2).squeeze()
l,m,n = data.shape
data = data.reshape(l,m,n/4,4).sum(axis=3).astype(np.float32)
nchan = n/4
print 'time spent quenching the frequency:', time.time() - now
now = time.time()

polX = []
polY = []


def runningmean(arr, windowsize):#running mean
    #init
    #fifo = deque(arr[:windowsize], maxlen = windowsize)
    #ave = sum(fifo) / windowsize
    ave = arr[:windowsize].mean()
    res = []
    arrlen = len(arr)
    for i in range(windowsize,  arrlen):
        a = arr[i]
        #head = fifo.pop()
        #fifo.append(a)
        head = arr[i-windowsize]
        res.append(head > ave)
        ave += (a - head)/windowsize
    for i in range(windowsize):
        a = arr[-1*i]
        #head = fifo.pop()
        #fifo.append(a)
        head = arr[arrlen + i - windowsize]
        res.append(head > ave)
        ave += (a - head)/windowsize
    return np.array(res).astype(np.uint8)


def func(i):
    dataS = data[:,:,i].flatten()
    rm = runningmean(dataS, WINDOWSIZE)
    #pb(i)
    return rm


#res = np.array(threadit(func, [[i] for i in range(nchan)])).T
onebitdata = np.array(threadit(func, [[i] for i in range(nchan)[::-1]])).T
print 'time spent medfilt the data:', time.time() - now

#from pylab import *
#imshow(onebitdata, aspect='auto')
#show()

#data = onebitdata.reshape(-1,4096,nchan)
#l,m,n = data.shape
#data = data.reshape(l*m,n)
#m,n = data.shape
#data = data.reshape(m/64, 64, n)
#data = data.mean(axis=1)

now = time.time()
data = np.packbits(onebitdata.reshape(-1,8)[:,::-1])
#data = onebitdata.flatten()
print 'bitpacking took:', time.time() - now
now = time.time()

from filwriter import filterbank_saver

fout = open('test.fil', 'wb')
ibeam = 0
nbeam = 1
nbits = 1
df *= -4
#print fout, sourcename, ibeam, nbeam, obsmjd, tstart, Azimuth, Zenith, ra, dec, nchan, nbits, tsamp, fch1, df
filterbank_saver(fout, sourcename, ibeam, nbeam, obsmjd, tstart, Azimuth, Zenith, src_raj, src_dej, nchan, nbits, tsamp, fch1, df)
fout.write(data)
fout.close()
print 'save file took:', time.time() - now

'''
from pylab import *
imshow(onebitdata.T, aspect='auto', origin='low')  #, extent=(0., float(duration), float(fmin), float(fmax)))
#xlabel('t (s)')
#ylabel('frequency (MHz)')
show()
'''


