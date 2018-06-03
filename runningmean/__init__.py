from ctypes import *
import os
import numpy as np
PATH, this_file = os.path.split(os.path.realpath(__file__))
lib = CDLL(PATH + "/runningmean.so")


def runningmean(data, windowsize):
    outdata = np.zeros(data.shape, dtype='u1')
    nsamp, nchan = data.shape
    #lib.run(data, nsamp, nchan, outdata, windowsize)
    lib.run(c_void_p(data.ctypes.data), c_int(nsamp), c_int(nchan), c_void_p(outdata.ctypes.data), c_int(windowsize))
    return outdata

