#!/usr/bin/python
"""
Script for plotting the Omega matrix.
"""

import theo_header, input_options, lib_file, error_handler
import numpy
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except:
    print "pylab/matplotlib not installed - plotting not possible"
    raise

class OmFrag_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def __init__(self, *args, **kwargs):
        self.state_list = []
        self.maxOm = 0.
        
        input_options.write_options.__init__(self, *args, **kwargs)

    ## \param fname file with the data produced in a previous analyze_tden.py run
    def read_OmFrag(self, fname='OmFrag.txt'):
        """
        Read the OmFrag.txt file written by analyze_tden.py
        """
        Ofile = open(fname, 'r')
        line = Ofile.next()
        
        numF = int(line)
        
        while True:
            try:
                line = Ofile.next()
            except StopIteration:
                break
            
            words = line.split()
            
            self.state_list.append({})
            state = self.state_list[-1]
            state['name'] = words[0]
            state['Om'] = words[1]
            state['OmFrag'] = numpy.zeros([numF, numF])
            
            for iel, el in enumerate(words[2:]):
                iF = iel % numF
                jF = iel / numF
                state['OmFrag'][iF, jF] = el
                
            self.maxOm = max(self.maxOm, state['OmFrag'].max())

    def OmFrag_input(self):
        plot_opts = ['Plot original values', 'Plot sqareroot scaled values']
        ichoice = self.ret_choose_list('Do you want to scale the values before plotting?', plot_opts, 1)
        self.write_option('plot_type', '%i'%ichoice)
        
        self.read_int('Resolution (dpi) for plotting', 'plot_dpi', 50)
        
        self.choose_list(
            'Coloring scheme for pseudocolor plots',
            'cmap',
        [
            ('Greys', 'Grey scale'),
            ('Oranges', 'Orange scale'),
            ('Blues', 'Blue scale'),
            ('RdBu', 'Red -> Blue'),
            ('YlGn', 'Yellow -> Green'),
            ('YlOrRd', 'Yellow -> Orange -> Red'),
            ('YlOrBr', 'Yellow -> Orange -> Brown'),
            ('RdYlBu', 'Red -> Yellow -> Blue')
        ], 'Greys'
        )
        
        self.read_int('Font size', 'fsize', 25)
        self.read_yn('Use the same scale for all plots', 'sscale', True)
        self.read_yn('Plot axes?', 'axis', False)
        self.read_yn('Plot colorbar for each individual plot?', 'cbar', False)
    
            
    def plot(self):
        hfname = 'OmFrag.html'
        hfile = lib_file.htmlfile(hfname)
        hfile.pre('Electron-hole correlation plots')
        hfile.write('<h2>Electron-hole correlation plots of the Omega matrices for the individual states.</h2>')
        
        htable = lib_file.htmltable(ncol=4)
        
        matplotlib.rc('font', size=self['fsize'])
        
        for state in self.state_list:
            if self['plot_type'] == 1:
                plot_arr = state['OmFrag']
            elif self['plot_type'] == 2:
                plot_arr = numpy.sqrt(state['OmFrag'])
            else:
                raise error_handler.ElseError('%i'%self['plot_type'], 'plot_type')
            
            if self['sscale']: vmax = self.maxOm
            else: vmax = state['OmFrag'].max()
            
            pylab.figure()
            pylab.pcolor(plot_arr, cmap=pylab.get_cmap(name=self['cmap']), vmin=0., vmax=vmax)
            
            if self['axis']: pylab.axis('on')                
            else: pylab.axis('off')
            
            if self['cbar']: pylab.colorbar()
            
            pname = 'pcolor_%s.png'%state['name']
            print "Writing %s ..."%pname
            pylab.savefig(pname, orientation='landscape',dpi=self['plot_dpi'])
            
            tel  = '<img src="%s", border="1" width="200">\n'%pname
            tel += '<br>%s'%state['name']
            htable.add_el(tel)
        
        if self['sscale']:
            pylab.figure()
            pylab.pcolor(numpy.zeros([1, 1]), cmap=pylab.get_cmap(name=self['cmap']), vmin=0., vmax=self.maxOm)
            pylab.axis('off')
            pylab.colorbar()
            pylab.savefig('cbar.png', dpi=self['plot_dpi'])
            
            tel  = '<img src="cbar.png", border="1" width="200">\n'
            tel += '<br>Scale'
            htable.add_el(tel)
            
        hfile.write(htable.ret_table())    
        hfile.post()
        
        print " HTML file %s containing the electron-hole correlation plots written."%hfname

def run_plot():
    Oopt = OmFrag_options('plot.in')        
    Oopt.read_OmFrag()
    
    Oopt.OmFrag_input()
    
    Oopt.plot()
    
    
if __name__ == '__main__':
    theo_header.print_header('Plot Omega matrices')

    run_plot()