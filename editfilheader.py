import numpy as np
import sys,os,glob
from struct import *
from ProgressBar import progressBar as PB

def update(filename):
    with open(filename, 'r+') as fin:
        header = fin.read(360)
        #print header
        i = header.find('tstart')
        #print int(header[i+5:i+13], 0)
        tstart = np.frombuffer(buffer(header, i+6, 8), dtype='float64')[0]

        wrongmjd = "%.5f" % tstart
        print 'tstart', wrongmjd
        wrongmjd = "%.5f" % (tstart / (1.+1./86400))
        wMJDi,wMJDs = wrongmjd.split('.')
        crrmjd = float(wMJDi) + float(wMJDs)/86400
        print 'mjd', crrmjd
        #mjdvalue = pack('d', crrmjd)
        #fin.seek(i+6)
        #fin.write(outvalue)
        #fin.close()

        #if foff == 0.4882812500:
            #newfoff = -0.4882812500
            #outvalue = pack('d', newfoff)
            #fin.seek(i+4)
            #fin.write(outvalue)
        #fin.close()

update(sys.argv[1])

#files = np.genfromtxt('wrongheaderfiles.txt', dtype="|S150")
#pb = PB(maxValue=len(files))
#for i,fil in enumerate(files[:]):
    ##files = glob.glob(subdir+'/*.fil')
    #update(fil)
    #pb(i)
