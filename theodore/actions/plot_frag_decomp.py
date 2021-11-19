#!/usr/bin/env python3
"""
Script for plotting fragment decomposition.
"""

from __future__ import print_function, division
from .actions import Action
import numpy
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    matplotlib = importer.lazy_import('matplotlib')
    pylab = importer.lazy_import('pylab')



class decomp_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def __init__(self, *args, **kwargs):
        self.state_list = []
        self.numF = 0 # number of fragments

        input_options.write_options.__init__(self, *args, **kwargs)
        self['colors'] = ['b', 'g', 'orange', 'red', 'gray']

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
        self.read_int('Font size', 'fsize', 8)
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

        pylab.figure(figsize=[len(self.state_list)+1, 5])

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
        pylab.xticks(ind, [state['name'] for state in self.state_list], rotation='vertical')

        pylab.subplots_adjust(bottom=0.2, left=0.2)

        fname = 'frag_decomp.%s'%self['output_format']
        pylab.savefig(fname, dpi=self['plot_dpi'])
        print(("File %s written."%fname))


class PlotFragDecomp(Action):
    name = 'plot_frag_decomp'

    _colt_description = 'Plot fragment decomposition of Omega matrix'

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
            '..error_handler': 'error_handler',
            'matplotlib': 'matplotlib',
            'pylab': 'pylab',
    })

    def run():
        matplotlib.use('Agg')
        theo_header.print_header(__class__._colt_description)
        opt = decomp_options('plot.in')
        opt.read_OmFrag()

        opt.decomp_input()

        opt.plot()
