#!/usr/bin/python
"""
Script for creating graphs from multiple directories, e.g. potential curves.
"""

import theo_header, input_options, error_handler
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except:
    print "pylab/matplotlib not installed - plotting not possible"
    raise

class plot_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def __init__(self, *args, **kwargs):
        input_options.write_options.__init__(self, *args, **kwargs)
        
    def plot_input(self):
        print "This script allows to combine information from several TheoDORE runs into one graph."
        print "   These jobs are assumed to be located in subdirectories of the current directory."
        
        sdirs = os.listdir('.')
        
        print "The following subdirectories were found:"
        self.print_list(sdirs)
        rstr = self.ret_str("Please enter the order in which they should appear, e.g. '1 2 4 3'")
        ilist = [int(idir) - 1 for idir in rstr.split()]
        dlist = [sdirs[idir] for idir in ilist]
        
        self.write_list('ana_dirs', dlist, "'%s'")
        
    def plot(self):
        pass
        
def run_plot():
    popt = plot_options('graph.in')
    
    popt.plot_input()
    
    popt.plot()
    
    popt.flush()
    
if __name__ == '__main__':
    theo_header.print_header('Graph plotting')
    
    run_plot()