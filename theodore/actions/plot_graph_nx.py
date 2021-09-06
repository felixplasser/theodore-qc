#!/usr/bin/env python3
"""
Script for creating graphs from multiple directories, e.g. potential curves.
"""

from __future__ import print_function, division
from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    lib_plot = importer.lazy_import_as('..lib_plot', 'lib_plot')
    lib_file = importer.lazy_import_as('..lib_file', 'lib_file')


class write_plot_options_nx(lib_plot.write_plot_options):
    def plot_input(self):
        self['ana_dirs']=[]
        self['ana_file']='nx.log'
        self['state_labels']=None
        self['dognu'] = False

        self.read_int('Number of states to plot', 'nstate', 10)
        self.read_float('Minimum time', 'tmin', 0.)
        self.read_float('Maximum time', 'tmax', 10000.)

        self.read_yn('Create plots using pylab?', 'doplots', True)
        if self['doplots']:
            self.read_int('Font size', 'fsize', 15)
            self.read_str("Format of output graphics files", "output_format", "png")

        self.read_yn('Print txt files with the information', 'dotxt', True)

    def read_data(self):
        self.data = []
        self.main_header = ''
        self.data_avail = []
        self.times = []
        self.act = []

        f = open(self['ana_file'])
        while True:
            try:
                line = next(f)
            except StopIteration:
                break

            if 'state       dE(eV)' in line:
                ddict = {}
                state_labels = []

                header = line.split()

                line = next(f) # ------
                line = next(f)

                for istate in range(self['nstate']):
                    if 'Finished' in line: break

                    words = line.split()
                    state_label = words[0]
                    ddict[state_label] = {}
                    pdict = ddict[state_label]
                    state_labels.append(state_label)

                    for i, prop in enumerate(header[1:]):
                        try:
                            pdict[prop] = float(words[i+1])
                        except ValueError:
                            pass
                    pdict['ddE'] = pdict['dE(eV)'] - ddict[state_labels[0]]['dE(eV)']
                    line = next(f)

            if 'FINISHING STEP' in line:
                words = line.split()
                time = float(words[4])

                if self['tmax'] <= time:
                    break
                elif self['tmin'] <= time:
                    self['ana_dirs'].append(time)
                    self.times.append(time)
                    act = int(words[-1]) - 2
                    ddict['act'] = ddict[state_labels[act]]
                    self.data.append(ddict)

        f.close()

        self.main_header = header + ['ddE']
        self['state_labels'] = state_labels + ['act']
        self['leg_labels'] = state_labels + ['act']


    def plot(self):
        """
        Create the plots.
        For this purpose, self.data has to be rearranged.
        """
        try:
            import matplotlib
            matplotlib.use('Agg')
            import pylab
        except:
            print("pylab/matplotlib not installed - plotting not possible")
            print("exiting ...")
            return

        hfname = 'graphs.html'
        hfile = lib_file.htmlfile(hfname)
        hfile.pre('Property graphs')
        htable = lib_file.htmltable(ncol=3)

        lfname = 'graphs.tex'
        lfile  = lib_file.latexfile(lfname)
        lfile.pre('Property graphs', graphicx=True)
        ltable = lib_file.latextable(ncol=2)

        #set1 = self.data[0][self['state_labels'][0]] # not used anywhere??

        matplotlib.rc('font', size=self['fsize'])

        for key in self.main_header[1:]:
            if key == 'fname': continue

            print('Plotting %s ...'%key)
            pylab.figure(figsize=(max(12,len(self.times)/100),6))

            for state in self['state_labels']:
                if state == 'act':
                    symb = 'ro'
                else:
                    symb = '-'

                ylist = []
                for iana_dir in range(len(self['ana_dirs'])):
                    try:
                        ylist.append(self.data[iana_dir][state][key])
                    except KeyError:
                        print(" ... not able to plot %s for %s."%(key, state))
                        break
                else:
                    pylab.plot(self.times, ylist, symb)

            pylab.title(key)
            pylab.ylabel(key)
            pylab.xlabel('time (fs)')

            pname = '%s.%s'%(key, self['output_format'])
            pylab.savefig(pname)

            tel  = '<img src="%s", border="1" width="%i">'%(pname, max(len(self.times)/2, 400))
            htable.add_el(tel)

            lel = "\\incplot{%s}"%pname
            ltable.add_el(lel)

        hfile.write(htable.ret_table())
        hfile.post(lvprt=1)

        lfile.write(ltable.ret_table())
        lfile.post(lvprt=1)


class PlotGraphNx(Action):

    name = 'plot_graph_nx'

    _colt_description = 'Graph plotting (Newton-X)'

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
            '..lib_plot': 'lib_plot',
            '..lib_file': 'lib_file',
    })


    def run():
        theo_header.print_header(title=__class__._colt_description)
        infilen = 'graph.in'

        popt = write_plot_options_nx(infilen)

        popt.plot_input()
        popt.read_data()

        if popt['doplots']: popt.plot()
        if popt['dotxt']:   popt.txt_files()
        if popt['dognu']:   popt.gnu_inp()
