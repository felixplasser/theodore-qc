#!/usr/bin/python
"""
Script for creating graphs from multiple directories, e.g. potential curves.
"""

import theo_header, input_options, error_handler, lib_file
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except:
    print "pylab/matplotlib not installed - plotting not possible"
    raise

class write_plot_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def __init__(self, *args, **kwargs):
        input_options.write_options.__init__(self, *args, **kwargs)
        
    def plot_input(self):
        print "This script allows to combine information from several TheoDORE runs into one graph."
        print "   These jobs are assumed to be located in subdirectories of the current directory."
        
        sdirs = [dirn for dirn in os.listdir('.') if os.path.isdir(dirn)]
        
        print "The following subdirectories were found:"
        self.print_list(sdirs)
        rstr = self.ret_str("Please enter the order in which they should appear, e.g. '1 2 4 3'")
        ilist = [int(idir) - 1 for idir in rstr.split()]
        dlist = [sdirs[idir] for idir in ilist]
        
        self.write_list('ana_dirs', dlist, "'%s'")
        
        self.read_str("Name of the file to analyze in each directory", "ana_file", "tden_summ.txt")
        
        rstr = self.ret_str("Labels of the states of interest as they appear in %s (separated by spaces)"%self['ana_file'])
        self.write_list('state_labels', rstr.split(), lformat="'%s'")
        
        self.read_int('Font size', 'fsize', 10)
        self.read_str("Format of output graphics files", "output_format", "png")
        
    def read_data(self):
        self.data = []
        self.main_header = ''
            
        for ana_dir in self['ana_dirs']:
            sfile = lib_file.summ_file(os.path.join(ana_dir, self['ana_file']))
            header = sfile.ret_header()
            ddict  = sfile.ret_ddict()
            
            if self.main_header == '': self.main_header = header
            self.data.append(ddict)            
            
        #print self.data
        
    def plot(self):
        hfname = 'graphs.html'
        hfile = lib_file.htmlfile(hfname)
        hfile.pre('Property graphs')
        
        htable = lib_file.htmltable(ncol=4)
        
        # TODO: read state labels
        set1 = self.data[0][self['state_labels'][0]]
        
        matplotlib.rc('font', size=self['fsize'])
        
        for key in self.main_header[1:]:
            if key == 'fname': continue
            
            print 'Plotting %s ...'%key
            pylab.figure(figsize=(6,4))
            
            for state in self['state_labels']:
                ylist = []
                for iana_dir in xrange(len(self['ana_dirs'])):
                    try:
                        ylist.append(self.data[iana_dir][state][key])
                    except KeyError:
                        print " ... not able to plot %s for %s."%(key, state)
                        break
                else:
                    pylab.plot(range(len(ylist)), ylist, 'x-', label=state)
            
            pylab.title(key)
            
            numx = len(self['ana_dirs'])
            pylab.xticks(xrange(numx), self['ana_dirs'], rotation=30)
            pylab.margins(0.20)
            pylab.subplots_adjust(bottom=0.15)
            pylab.xlim((-0.5, numx+1.5))
            
            pylab.ylabel(key)
            pylab.legend()
            
            pname = '%s.%s'%(key, self['output_format'])
            pylab.savefig(pname)
            
            tel  = '<img src="%s", border="1" width="400">'%pname
            htable.add_el(tel)
                
        hfile.write(htable.ret_table()) 
        hfile.post()
        
        print " HTML file %s containing the property graphs written."%hfname

class read_plot_options(input_options.read_options):
    def set_defaults(self):
        self['ana_dirs']=None
        self['ana_file']='tden_summ.txt'
        self['state_labels']=None
        self['fsize']=10
        self['output_format']='png'
        
def run_plot():
    infilen = 'graph.in'
    
    popt = write_plot_options(infilen)
    ropt = read_plot_options(infilen, False)

    if ropt.init == 0:
        copy = popt.ret_yn('Found %s. Use this file directly rather than performing an interactive input?'%infilen, False)
    else:
        copy = False    
    
    if copy:
        popt.copy(ropt)
    else:        
        popt.plot_input()
    
    popt.read_data()
    
    popt.plot()
    
    if not copy:
        popt.flush()
    
if __name__ == '__main__':
    theo_header.print_header('Graph plotting')
    
    run_plot()