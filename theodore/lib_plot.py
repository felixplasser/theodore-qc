"""
Library file with some routines for plotting graphs.
"""

from __future__ import print_function, division

from . import input_options, error_handler, lib_file
import os

class write_plot_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def __init__(self, *args, **kwargs):
        input_options.write_options.__init__(self, *args, **kwargs)
        self['symb'] = 'x-'

    def plot_input(self):
        """
        Read input from command line.
        """
        print("This script allows to combine information from several TheoDORE runs into one graph.")
        print("   These jobs are assumed to be located in subdirectories of the current directory.")

        sdirs = sorted([dirn for dirn in os.listdir('.') if os.path.isdir(dirn)])

        print("The following subdirectories were found:")
        self.print_list(sdirs)
        rstr = self.ret_str("Please enter the order in which they should appear, e.g. '1 2 4 3'")
        ilist = [int(idir) - 1 for idir in rstr.split()]
        dlist = [sdirs[idir] for idir in ilist]

        self.write_list('ana_dirs', dlist, "'%s'")

        self.read_str("Name of the file to analyze in each directory", "ana_file", "tden_summ.txt")

        labstr = self.ret_str("Labels of the states of interest as they appear in %s (separated by spaces)"%self['ana_file'])
        self.write_list('state_labels', labstr.split())

        rstr = self.ret_str("State labels for the legend", labstr)
        self.write_list('leg_labels', rstr.split())

        self.read_yn('Create plots using pylab?', 'doplots', True)
        if self['doplots']:
            self.read_int('Font size', 'fsize', 10)
            self.read_str("Format of output graphics files", "output_format", "png")

        self.read_yn('Create input for gnuplot?', 'dognu', False)
        if self['dognu']:
            self['dotxt'] = True
        else:
            self.read_yn('Print txt files with the information', 'dotxt', True)

    def read_data(self):
        """
        Read the data from the individual directories into individual dictionaries.
        Arranged as:
        [ana_dir] - {state_label} - {key}
        """
        self.data = []
        self.main_header = ''
        self.data_avail = []

        for ana_dir in self['ana_dirs']:
            sfile = lib_file.summ_file(os.path.join(ana_dir, self['ana_file']))
            header = sfile.ret_header()
            ddict  = sfile.ret_ddict()

            if self.main_header == '': self.main_header = header
            self.data.append(ddict)

        #print self.data

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
            raise

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
            pylab.figure(figsize=(6,4))

            for i, state in enumerate(self['state_labels']):
                ylist = []
                for iana_dir in range(len(self['ana_dirs'])):
                    try:
                        ylist.append(self.data[iana_dir][state][key])
                    except KeyError:
                        print(" ... not able to plot %s for %s."%(key, state))
                        break
                else:
                    pylab.plot(list(range(len(ylist))), ylist, self['symb'], label=self['leg_labels'][i])

            pylab.title(key)
            pylab.ylabel(key)

            numx = len(self['ana_dirs'])
            pylab.xticks(range(numx), self['ana_dirs'], rotation=30)
            #pylab.margins(0.20)
            pylab.subplots_adjust(bottom=0.15)
            pylab.xlim((-0.5, numx+1.5))
            pylab.legend()

            pname = '%s.%s'%(key, self['output_format'])
            pylab.savefig(pname)

            tel  = '<img src="%s", border="1" width="400">'%pname
            htable.add_el(tel)

            lel = "\\incplot{%s}"%pname
            ltable.add_el(lel)

        hfile.write(htable.ret_table())
        hfile.post(lvprt=1)

        lfile.write(ltable.ret_table())
        lfile.post(lvprt=1)

    def txt_files(self):
        """
        Create compact text files that contain all the required info.
        """
        for key in self.main_header[1:]:
            if key == 'fname': continue

            found_data = False

            fname = '%s.txt'%key
            print('Writing %s ...'%fname)

            wf = open(fname, 'w')

            wf.write('%10s'%'x')
            for state in self['leg_labels']:
                wf.write('%10s'%state)

            wf.write('\n')
            for idir, ana_dir in enumerate(self['ana_dirs']):
                wf.write('%10s'%ana_dir)

                for state in self['state_labels']:
                    try:
                        wf.write('%10.5f'%self.data[idir][state][key])
                        found_data = True
                    except KeyError:
                        #wf.write('    -     ')
                        wf.write('    nan   ')

                wf.write('\n')

            wf.close()

            if found_data: self.data_avail.append(key)

    def gnu_inp(self):
        """
        Create input for gnuplot.
        """
        #TODO: One could take all the settings out of the loop and then write the plot commands explicitly.

        propstr = ""
        for key in self.main_header[1:]:
            if key in self.data_avail:
                propstr += " " + key

        ticstr = "set xtics ("
        for i, ad in enumerate(self['ana_dirs']):
            if i > 0: ticstr += ", "
            ticstr += "'%s' %i"%(ad, i+1)
        ticstr += ")\n"

        plotstr = "plot"
        for i, leg in enumerate(self['leg_labels']):
            if i > 0: plotstr += ","
            plotstr += ' name.".txt" using %i title "%s"'%(i+2, leg)
        plotstr += "\n"

        fname = 'graph.gnuplot'
        with open(fname, 'w') as wf:
            wf.write("set terminal pngcairo size 800,600 enhanced font 'Helvetica,15' linewidth 3\n\n")
            wf.write('set xlabel "Directory"\n')
            wf.write('set xrange [%.1f:%.1f]\n'%(0.5, len(self['ana_dirs']) + 0.5))
            wf.write(ticstr)
            wf.write('set xtics nomirror rotate by -45\n')
            wf.write('set style data linespoints\n\n')
            wf.write('do for [name in "%s"] {\n'%propstr)
            wf.write('  set output name.".png"\n\n')
            wf.write('set ylabel name\n')
            wf.write(plotstr)
            wf.write('}\n')
            wf.write('quit\n')

        print(("gnuplot script %s created."%fname))

class read_plot_options(input_options.read_options):
    def set_defaults(self):
        self['ana_dirs']=None
        self['ana_file']='tden_summ.txt'
        self['state_labels']=None
        self['leg_labels']=None
        self['fsize']=10
        self['output_format']='png'
        self['doplots']=True
        self['dotxt']=True
        self['dognu']=True
