#!/usr/bin/env python2
"""
Script for plotting the Omega matrix as a bar plot.
This creates a latex/pgfplots input file.
Author: Felix Plasser, Sebastian Mai
"""

import theo_header, input_options, lib_file, error_handler
import numpy
import os

class Om_bar_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def read_OmFrag(self, fname='OmFrag.txt'):
        """
        Read the OmFrag.txt file written by analyze_tden.py
        """
        self.state_list = []
        Ofile = open(fname, 'r')
        line = Ofile.next()

        self.numF = int(line)

        while True:
            try:
                line = Ofile.next()
            except StopIteration:
                break

            words = line.split()

            self.state_list.append({'name':words[0],
                'Om':words[1],
                'OmFrag':map(float, words[2:])})
        self.numSt = len(self.state_list)

    def Om_bar_input(self):
        self.read_str("Name of the file with the tden information", "tdenfile", "tden_summ.txt", autocomp=True)

        if self.numSt < 20:
            defw = 7.
        else:
            defw = 15.
        self.read_float("Width of the plot (cm)", 'width', defw)

        self.comps = []
        print "Please enter the different excitation components to be plotted"
        print "    - leave empty to finish"
        for icomp in range(1,1000):
            rstr = self.ret_str('Name of component %i'%icomp)
            if rstr == '':
                print " ... component input finished."
                break

            color = self.ret_str('Color for plotting')

            print "\n *** Fragment pairs belonging to %s ***"%rstr
            print "  Enter two indices between 1 and %i, separated by spaces"%self.numF
            print "  Leave empty to finish"
            cols = [] # Columns in the OmFrag.txt file
            for jpair in range(1,1000):
                ehstr = self.ret_str("Hole/electron indices for pair %i"%jpair)
                if ehstr == '':
                    print " ... switching to next component."
                    break

                words = ehstr.split()
                jh = int(words[0]) - 1
                je = int(words[1]) - 1
                cols.append(jh * self.numF + je)

            self.comps.append({'name':rstr, 'color':color, 'cols':cols})

    def Om_bar_data(self):
        """
        Create file with the data to be plotted.
        """
        sf = lib_file.summ_file(self['tdenfile'])
        f = open('Om_bar_data.txt', 'w')
        f.write("Header\n")

        ist = 1
        for state in self.state_list:
            dE = sf.ddict[state['name']]['dE(eV)']
            try:
                osc = sf.ddict[state['name']]['f']
            except KeyError:
                osc = 0.
            f.write("%3i %10s % .5f % .5f"%(ist, state['name'], dE, osc))
            ist += 1
            for comp in self.comps:
                wt = sum(state['OmFrag'][icol] for icol in comp['cols'])
                f.write(" %.6f"%wt)
            f.write("\n")
        f.close()

    def make_tex(self):
        lfile = lib_file.latexfile('Om_bars.tex')
        lfile.write(self.pre())

        # Bar graphs with state characters
        lfile.write("""\\begin{axis}[
      anchor=north west, at={(0,0)},
      xlabel={State}, ylabel={Character},
      ytick={0.0,0.2,...,0.8},\n""")
        lfile.write("ymin=-0.0, ymax=1,\n")
        lfile.write("ybar stacked, bar width=%.3f cm]\n"%(self['width']/self.numSt/1.2))

        for icomp, comp in enumerate(self.comps):
            lfile.write("\\addplot[ybar, draw=none, fill=%s] table[x index=0, y index=%i] {Om_bar_data.txt};\n"%(comp['color'], icomp+4))
            lfile.write("\\addlegendentry{%s};\n"%comp['name'])
        lfile.write('\end{axis}\n\n')

        # Graphs with energies
        lfile.write("""\\begin{axis}[
      anchor=south west, at={(0,0)},
      height=2.0cm,
      xlabel={}, ylabel={Energy (eV)}
      %, colorbar % uncomment for colorbar as legend
      ]\n""")

        lfile.write("\\addplot[mark=-, only marks, mark options={scale=2, line width=2pt}]  table[x index=0, y index=2] {Om_bar_data.txt};\n\n")
        lfile.write("% Uncomment to plot oscillator strengths as shading\n")
        lfile.write("%\\addplot+[solid, mesh, point meta=explicit, no markers, line width=2cm, shader=flat corner] table[x expr=\\thisrowno{0}-0.5, y expr=5, meta expr=\\thisrowno{3}] {Om_bar_data.txt};\n")
        lfile.write('\end{axis}\n')

        # Finish
        lfile.write('\end{tikzpicture}\n')
        lfile.post(lvprt=1)
        print "  -> Create plots using: pdflatex %s"%lfile.name

    def pre(self):
        """
        Code for beginning of tex file.
        """
        str = """\documentclass{standalone}
\usepackage{xcolor}
\usepackage{tikz, pgfplots}
\usetikzlibrary{calc, decorations, plotmarks, fit, positioning, shapes.geometric}
\pgfplotsset{compat=1.4}

% ===========================================================================

% settings for the coordinate system
\pgfplotsset{
every linear axis/.append style={
"""
        str += 'width=%.3f cm, height=2.0cm,'%self['width']
        str += "xmin=0, xmax=%i,\n"%(self.numSt+1)
        str += """xtick={-10}, xticklabels={},
  scale only axis, axis on top,
  yticklabel style={
    /pgf/number format/fixed,
    /pgf/number format/zerofill,
    /pgf/number format/precision=1},
  legend style={at={(1.01,0.00)}, anchor=south west, draw=none, fill=black!5, inner sep=0.5pt, outer sep=0.5pt},
  legend columns=1, legend cell align=left},
}
\pgfplotsset{colormap={CI}{color=(white); color=(green!50); color=(blue!70);}}

% ===========================================================================

\\begin{document}
\\begin{tikzpicture}
"""
        return str

def run_plot():
    opt = Om_bar_options('plot.in')
    opt.read_str("Name of the file with the Omega matrix entries", "OmFfile", "OmFrag.txt", autocomp=True)
    opt.read_OmFrag(opt['OmFfile'])
    opt.Om_bar_input()
    opt.Om_bar_data()
    opt.make_tex()


if __name__ == '__main__':
    theo_header.print_header('Plot Omega matrices as bar graphs')

    run_plot()
