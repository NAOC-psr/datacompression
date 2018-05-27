from ctypes import *
import os
PATH, this_file = os.path.split(os.path.realpath(__file__))
lib = CDLL(PATH + "/runningmean.so")


def runningmean(data, outdata, windowsize):
    nsamp, npol, nchan = data.shape
    #lib.run(data, nsamp, nchan, outdata, windowsize)
    lib.run(c_void_p(data.ctypes.data), c_int(nsamp), c_int(nchan), c_void_p(outdata.ctypes.data), c_int(windowsize))

