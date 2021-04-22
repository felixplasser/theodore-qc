#!/usr/bin/env python3
"""
Script for plotting the Omega matrix as a pseudocolor matrix plot.
"""

from __future__ import print_function, division
from .. import theo_header, input_options, lib_file, error_handler
import numpy
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except:
    print("pylab/matplotlib not installed - plotting not possible")
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
        line = next(Ofile)

        numF = int(line)

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
            state['OmFrag'] = numpy.zeros([numF, numF])

            for iel, el in enumerate(words[2:]):
                iF = iel %  numF
                jF = iel // numF
                state['OmFrag'][iF, jF] = el

            self.maxOm = max(self.maxOm, state['OmFrag'].max())

    def OmFrag_input(self):
        plot_opts = ['Plot original values', 'Plot sqareroot scaled values']
        ichoice = self.ret_choose_list('Do you want to scale the values before plotting?', plot_opts, 1)
        self.write_option('plot_type', ichoice)

        self.read_int('Resolution (dpi) for plotting', 'plot_dpi', 200)

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
            ('RdYlBu', 'Red -> Yellow -> Blue'),
            ('binary', 'White -> Black'),
            ('bone_r', 'White -> Blue -> Black'),
            ('pink_r', 'White -> Pink -> Black'),
            ('hot_r',  'White -> Yellow -> Red -> Black')
        ], 'Greys'
        )

        self.read_int('Font size', 'fsize', 10)
        self.read_str("Format of output graphics files", "output_format", "png", autocomp=False)
        self.read_yn('Use the same scale for all plots', 'sscale', True)
        if self['sscale']:
            self.read_float('Minimal value to plot', 'vmin', 0.)
            self.read_float('Maximal value to plot', 'vmax', self.maxOm)
        self.read_yn('Plot frame?', 'axis', True)
        if self['axis']:
            self.read_yn('Axis with tick labels?', 'ticks', False)
            if self['ticks']:
                self.read_yn('Use custom tick labels?', 'cticks', False)
                if self['cticks']:
                    self['xticks'] = self.ret_str('Enter x-tick labels (separated by spaces)').split()
                    self['yticks'] = self.ret_str('Enter y-tick labels (separated by spaces)').split()
        self.read_yn('Draw grid?', 'grid', True)
        self.read_yn('Plot colorbar for each individual plot?', 'cbar', False)


    def plot(self):
        hfname = 'OmFrag.html'
        hfile = lib_file.htmlfile(hfname)
        hfile.pre('Electron-hole correlation plots')
        hfile.write('<h2>Electron-hole correlation plots of the Omega matrices for the individual states.</h2>')

        htable = lib_file.htmltable(ncol=4)

        matplotlib.rc('font', size=self['fsize'])

        if self['grid']:
            edgecolors='k'
        else:
            edgecolors=None

        for state in self.state_list:
            if self['plot_type'] == 1:
                plot_arr = state['OmFrag']
            elif self['plot_type'] == 2:
                plot_arr = numpy.sqrt(state['OmFrag'])
            else:
                raise error_handler.ElseError(str(self['plot_type']), 'plot_type')

            if self['sscale']:
                vmin = self['vmin']
                vmax = self['vmax']
            else:
                vmin = 0.
                vmax = state['OmFrag'].max()

            # Completely delete the small elements
            # for x in numpy.nditer(plot_arr, op_flags = ['readwrite']):
            #     if x < vmin:
            #         x[...] = -1. # numpy.nan

            pylab.figure(figsize=(2,2))
            pylab.pcolor(plot_arr, cmap=pylab.get_cmap(name=self['cmap']), vmin=vmin, vmax=vmax, edgecolors=edgecolors)

            # *** Different colouring of different parts ***
            # frag_lists = [[0, 2, 4, 6], [1, 3, 5]]
            # cmaps = ['Reds', 'Blues']
            # OmDim = len(plot_arr)
            # for frag in frag_lists:
            #     tmp_arr = numpy.array([[numpy.nan for i in range(OmDim)] for j in range(OmDim)])
            #     for i in frag:
            #         tmp_arr[i,i] = plot_arr[i,i]
            #     pylab.pcolor(tmp_arr, cmap=pylab.get_cmap(cmaps.pop(0)), vmin=0., vmax=vmax, edgecolors=edgecolors)

            if self['axis']:
                pylab.axis('on')
                if self['ticks']:
                    pylab.tick_params(which='both', length=0)
                    if self['cticks']:
                        pylab.xticks([x + 0.5 for x in range(len(plot_arr))], self['xticks'])
                        pylab.yticks([y + 0.5 for y in range(len(plot_arr))], self['yticks'])
                    else:
                        pylab.xticks([x + 0.5 for x in range(len(plot_arr))], [x + 1 for x in range(len(plot_arr))])
                        pylab.yticks([y + 0.5 for y in range(len(plot_arr))], [y + 1 for y in range(len(plot_arr))])
                else:
                    pylab.xticks([])
                    pylab.yticks([])
            else:
                pylab.axis('off')

            if self['cbar']: pylab.colorbar()

            pname = 'pcolor_%s.%s'%(state['name'], self['output_format'])
            print("Writing %s ..."%pname)
            pylab.tight_layout()
            pylab.savefig(pname, dpi=self['plot_dpi'])
            pylab.close()

            tel  = '<img src="%s", border="1" width="200">\n'%pname
            tel += '<br>%s'%state['name']
            htable.add_el(tel)

        # create a plot with the e/h axes and optionally the scale
        pylab.figure(figsize=(3,2))
        matplotlib.rc('font', size=14)
        ax = pylab.axes()
        ax.arrow(0.15, 0.15, 0.5, 0., head_width=0.05, head_length=0.1, fc='r', ec='r')
        ax.text(0.20, 0.03, 'hole', color='r')
        ax.arrow(0.15, 0.15, 0., 0.5, head_width=0.05, head_length=0.1, fc='b', ec='b')
        ax.text(0.02, 0.20, 'electron', rotation='vertical', color='b')

        pylab.axis('off')
        if self['sscale']:
            pylab.savefig('axes_no.%s'%self['output_format'], dpi=self['plot_dpi'])
#            pylab.figure(figsize=(2,2))

            pylab.pcolor(numpy.zeros([1, 1]), cmap=pylab.get_cmap(name=self['cmap']), vmin=self['vmin'], vmax=self['vmax'])

            pylab.colorbar()

        pylab.savefig('axes.%s'%self['output_format'], dpi=self['plot_dpi'])

        tel  = '<img src="axes.%s", border="1" width="200">\n'%self['output_format']
        tel += '<br>Axes / Scale'
        htable.add_el(tel)

        hfile.write(htable.ret_table())
        hfile.post()

        print(" HTML file %s containing the electron-hole correlation plots written."%hfname)

def run_plot():
    Oopt = OmFrag_options('plot.in')
    Oopt.read_OmFrag()

    Oopt.OmFrag_input()

    Oopt.plot()


def plot_omfrag():
    theo_header.print_header('Plot Omega matrices')
    run_plot()
