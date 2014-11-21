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
        
        self.read_str("Name of the file to analyze in each directory", "ana_file", "tden_summ.txt")
        
        rstr = self.ret_str("Labels of the states of interest as they appear in %s (separated by spaces)"%self['ana_file'])
        self.write_list('state_labels', rstr.split(), lformat="'%s'")
        
    def read_data(self):
        self.data = []
        self.main_header = ''
            
        for ana_dir in self['ana_dirs']:
            self.data.append({})
            f = open(os.path.join(ana_dir, self['ana_file']), 'r')
            
            header = f.next().split()
            if self.main_header == '': self.main_header = header
            f.next()
            
            while True:
                try:
                    line = f.next()
                except StopIteration:
                    break
                
                words = line.split()
                self.data[-1][words[0]] = {}
                pdict = self.data[-1][words[0]]
                
                for i, prop in enumerate(header[1:]):
                    try:
                        pdict[prop] = float(words[i+1])
                    except ValueError:
                        pass
            
        #print self.data
        
    def plot(self):
        set1 = self.data[0][self['state_labels'][0]]
        
        for key in self.main_header[1:]:
            print 'Plotting %s ...'%key
            pylab.figure(figsize=(6,4))
            
            for state in self['state_labels']:
                ylist = []
                for iana_dir in xrange(len(self['ana_dirs'])):
                    print self.data[iana_dir][state]
                    #ylist.append(self.data[iana_dir][state])
                    
                raise error_handler.NIError()
                
        
def run_plot():
    popt = plot_options('graph.in')
    
    popt.plot_input()
    
    popt.read_data()
    
    popt.plot()
    
    popt.flush()
    
if __name__ == '__main__':
    theo_header.print_header('Graph plotting')
    
    run_plot()