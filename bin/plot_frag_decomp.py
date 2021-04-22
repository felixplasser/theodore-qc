#!/usr/bin/env python3
"""
Script for plotting fragment decomposition.
"""

from __future__ import print_function, division
from .. import theo_header, input_options, error_handler
import numpy

try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except:
    print("pylab/matplotlib not installed - plotting not possible")
    raise

class decomp_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def __init__(self, *args, **kwargs):
        self.state_list = []
        self.numF = 0 # number of fragments

        input_options.write_options.__init__(self, *args, **kwargs)
        self['colors'] = ['b', 'g', 'orange', 'red']

    ## \param fname file with the data produced in a previous analyze_tden.py run
    def read_OmFrag(self, fname='OmFrag.txt'):
        """
        Read the OmFrag.txt file written by analyze_tden.py
        """
        Ofile = open(fname, 'r')
        line = next(Ofile)

        self.numF = int(line)

        while True:
            try:
                line = next(Ofile)
            except StopIteration:
                break

            words = line.split()

            self.state_list.append({})
            state = self.state_list[-1]
            state['name'] = words[0]
            state['Om'] = words[1]
            state['OmFrag'] = numpy.zeros([self.numF, self.numF], float)

            for iel, el in enumerate(words[2:]):
                iF = iel %  self.numF
                jF = iel // self.numF
                state['OmFrag'][iF, jF] = el # max(float(el), 0.)

    def decomp_input(self):

        self.read_float('Relative width of the bars', 'barwidth', 0.75)
        self.read_int('Resolution (dpi) for plotting', 'plot_dpi', 200)
        self.read_int('Font size', 'fsize', 5)
        self.read_str("Format of output graphics files", "output_format", "png", autocomp=False)
        self.labels = []
        for iF in range(self.numF):
            self.labels.append(self.ret_str("Label for fragment %i"%(iF+1), 'F%i'%(iF+1)))

    def plot(self):
        matplotlib.rc('font', size=self['fsize'])

        hpops = numpy.zeros([len(self.state_list),self.numF])
        epops = numpy.zeros([len(self.state_list),self.numF])
        ind = numpy.arange(len(self.state_list))

        for istate, state in enumerate(self.state_list):
            hpops[istate, :] = -numpy.sum(state['OmFrag'], 0)
            epops[istate, :] =  numpy.sum(state['OmFrag'], 1)

        #print hpops
        #print epops

        pylab.figure(figsize=[0.5 * len(self.state_list)+1,3])

        barkwargs = {'width':self['barwidth']}

        ebottom = numpy.zeros(len(self.state_list))
        hbottom = numpy.zeros(len(self.state_list))
        for iF in range(self.numF):
            pylab.bar(ind, epops[:, iF], bottom=ebottom, color=self['colors'][iF], label=self.labels[iF], **barkwargs)
            ebottom += epops[:, iF]

            pylab.bar(ind, hpops[:, iF], bottom=hbottom, color=self['colors'][iF], **barkwargs)
            hbottom += hpops[:, iF]

        pylab.xlabel('Excited states')
        pylab.ylabel('  Hole   <---> Electron')
        pylab.plot([ind[0]-0.5, ind[-1]+0.5+self['barwidth']], [0.,0.], 'k-')

        pylab.legend()
        pylab.xticks([])

        fname = 'frag_decomp.%s'%self['output_format']
        pylab.savefig(fname, dpi=self['plot_dpi'])
        print(("File %s written."%fname))

def run_plot():
    opt = decomp_options('plot.in')
    opt.read_OmFrag()

    opt.decomp_input()

    opt.plot()


def plot_frag_decomp():
    theo_header.print_header('Plot fragment decomposition')
    run_plot()
