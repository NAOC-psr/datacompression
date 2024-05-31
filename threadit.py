"""
A module for threading computation on a multi-cpu machine

by Weiwei Zhu
June 2013
"""
import sys
import traceback
import multiprocessing as MP
from inspect import getargspec
num_cpus = max(1, MP.cpu_count() - 1)


def threadit(func, arglist, OnOffSwitch={'state':False}, num_threads=50):
    """
    A wrapper for multi-threading any function (func) given a argument list (arglist). The OnOffSwitch is a flag that got set to True when a progress is already running in a thread. It would not spam more threads when the flag is set to True.
    """
    from threading import Thread
    resultdict = {}
    class thread(Thread):
        def __init__(self, threadID, func, args):
            self.threadID = threadID
            self.func = func
            self.args = args
            Thread.__init__(self)
        def run(self):
            res = self.func(*self.args)
            resultdict.update({self.threadID : res})

    if OnOffSwitch['state'] == False or len(arglist) <=3:
        OnOffSwitch['state'] = True
        for i,arg in enumerate(arglist):
            t = thread(i, func, arg)
            t.start()
            t.join()
        OnOffSwitch['state'] = False
    else:
        resultdict = {}
        for i in range(len(arglist)):
            resultdict.update({i:func(*(arglist[i]))})
    return [ resultdict[i] for i in range(len(resultdict))]
    
def spam(func):
    args = getargspec(func)[0]
    print('func args:', args)
    if not 'worker_id' in args:
        print('without worker_id')
        kws = {'worker_id':0}
        def wrap_func(*args, **kws):
            return func(*args)
        return wrap_func
    else:
        print('with worker_id')
        return func


def spamit(func, arglist, OnOffSwitch={'state':False}, num_threads=50):
    """
    A wrapper for multi-processing any function (func) given a argument list (arglist). The OnOffSwitch is a flag that got set to True when a progress is already running in a thread. It would not spam more threads when the flag is set to True.
    """
    num_workers = min(num_threads, num_cpus)
    def worker(worker_id, q,retq, pipe, func, arglist):
        while True:
            idx = q.get()
            if idx is not None:
                try:
                    retq.put({idx:func(*(arglist[idx]))})
                except:
                    except_type, except_class, tb = sys.exc_info()
                    pipe.send((except_type, except_class, traceback.extract_tb(tb)))

                    retq.put(None)
            else:
                break
            q.task_done()
        q.task_done()
    #print func.__name__, ' OnOffSwitch:', OnOffSwitch['state']
    if OnOffSwitch['state'] == False or len(arglist) <=3:
        #if no threading is already running or the number of jobs to spaw is smaller than 3, don't thread it.
        OnOffSwitch['state'] = True
        q=MP.JoinableQueue()
        to_child, to_self = MP.Pipe()
        retq=MP.Queue()
        procs = []
        for i in range(num_workers):
            p = MP.Process(target=worker, args=(i, q, retq, to_self, func, arglist))
            p.daemon = True
            p.start()
            procs.append(p)

        for i in range(len(arglist)):
            q.put(i)

        for p in range(num_workers):
            q.put(None)
        q.join()
        resultdict = {}
        for i in range(len(arglist)):
            res = retq.get()
            if not res == None:
                resultdict.update(res)
            else:
                exc_info = to_child.recv()
                print(exc_info)
                #raise exc_info[1]
        for p in procs:
            p.join()
        OnOffSwitch['state'] = False
        #return resultdict
    else:
        resultdict = {}
        for i in range(len(arglist)):
            resultdict.update({i:func(*(arglist[i]))})
        #return resultdict
    return [ resultdict[i] for i in range(len(resultdict))]


def workerpool(func, arglist, OnOffSwitch={'state':False}, num_threads=50):
    """
    A wrapper for multi-processing any function (func) given a argument list (arglist). The OnOffSwitch is a flag that got set to True when a progress is already running in a thread. It would not spam more threads when the flag is set to True.
    """
    num_workers = min(num_threads, num_cpus)
    def worker(worker_id, q,retq, pipe, func, arglist):
        while True:
            idx = q.get()
            if idx is not None:
                try:
                    retq.put({idx:func(*(arglist[idx]), worker_id=worker_id)})
                except:
                    except_type, except_class, tb = sys.exc_info()
                    pipe.send((except_type, except_class, traceback.extract_tb(tb)))

                    retq.put(None)
            else:
                break
            q.task_done()
        q.task_done()
    #print func.__name__, ' OnOffSwitch:', OnOffSwitch['state']
    if OnOffSwitch['state'] == False or len(arglist) <=3:
        #if no threading is already running or the number of jobs to spaw is smaller than 3, don't thread it.
        OnOffSwitch['state'] = True
        q=MP.JoinableQueue()
        to_child, to_self = MP.Pipe()
        retq=MP.Queue()
        procs = []
        for i in range(num_workers):
            p = MP.Process(target=worker, args=(i, q, retq, to_self, func, arglist))
            p.daemon = True
            p.start()
            procs.append(p)

        for i in range(len(arglist)):
            q.put(i)

        for p in range(num_workers):
            q.put(None)
        q.join()
        resultdict = {}
        for i in range(len(arglist)):
            res = retq.get()
            if not res == None:
                resultdict.update(res)
            else:
                exc_info = to_child.recv()
                print(exc_info)
                #raise exc_info[1]
        for p in procs:
            p.join()
        OnOffSwitch['state'] = False
        #return resultdict
    else:
        resultdict = {}
        for i in range(len(arglist)):
            resultdict.update({i:func(*(arglist[i]))})
        #return resultdict
    return [ resultdict[i] for i in range(len(resultdict))]




class hamsterball(object):
    
    def __init__(self):
        self.cwd = os.getcwd()
        self.tmpdir = tempfile.mkdtemp(prefix=self.cwd+'/tmp')

    def run(self, infiles, outfiles, proc, args):
        os.chdir(self.tmpdir)
        for infile in infiles: 
            os.symlink(self.cwd+'/'+infile, self.tmpdir+'/'+infile)
            #os.system("ln -s %s %s" % (self.cwd+'/'+infile, self.tmpdir+'/'+infile))
            #print self.cwd+'/'+infile, self.tmpdir+'/'+infile
        result = proc(*args)
        for infile in infiles:
            os.remove(self.tmpdir+'/'+infile)
        for outpat in outfiles:
            for outfile in glob.glob(outpat): 
                shutil.move(outfile, self.cwd+'/'+outfile)
        os.chdir(self.cwd)
        self.cleanup()
        return result

    def cleanup(self):
        shutil.rmtree(self.tmpdir)

def hamster(infiles, outfiles, proc, args):
    newball = hamsterball()
    return newball.run(infiles, outfiles,  proc, args)




class hamsterball(object):
    
    def __init__(self):
        self.cwd = os.getcwd()
        self.tmpdir = tempfile.mkdtemp(prefix=self.cwd+'/tmp')

    def run(self, infiles, outfiles, proc, args):
        os.chdir(self.tmpdir)
        for infile in infiles: 
            os.symlink(self.cwd+'/'+infile, self.tmpdir+'/'+infile)
            #os.system("ln -s %s %s" % (self.cwd+'/'+infile, self.tmpdir+'/'+infile))
            #print self.cwd+'/'+infile, self.tmpdir+'/'+infile
        result = proc(*args)
        for infile in infiles:
            os.remove(self.tmpdir+'/'+infile)
        for outpat in outfiles:
            for outfile in glob.glob(outpat): 
                shutil.move(outfile, self.cwd+'/'+outfile)
        os.chdir(self.cwd)
        self.cleanup()
        return result

    def cleanup(self):
        shutil.rmtree(self.tmpdir)

def hamster(infiles, outfiles, proc, args):
    newball = hamsterball()
    return newball.run(infiles, outfiles,  proc, args)

