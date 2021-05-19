import os
import argparse
# Python 2/3 compatibility
try:
    from time import process_time
except ImportError:
    from time import clock as process_time # clock is depricated and will be removed in 3.8
try:
    from time import perf_counter
except ImportError:
    from time import time as perf_counter



def isfile(filename):
    if not os.path.exists(filename):
        raise SystemExit('Input file %s not found!\nPlease create this file using theoinp or specify its location using -ifile\n' % filename)
    return filename


def get_ifile_commandline():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('-ifile', default='dens_ana.in', help='name of the input file')
    args = parser.parse_args()
    return isfile(args.ifile)


def timeit(func):
    def _wrapper(*args, **kwargs):
        (tc, tt) = (process_time(), perf_counter())
        result = func(*args, **kwargs)
        print("CPU time: % .1f s, wall time: %.1f s"%(process_time() - tc, perf_counter() - tt))
        return result
    return _wrapper
