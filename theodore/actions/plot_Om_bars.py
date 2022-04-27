#!/usr/bin/env python3
"""
Script for plotting the Omega matrix as a bar plot.
This creates a latex/pgfplots input file.
Author: Felix Plasser, Sebastian Mai
"""

from __future__ import print_function, division
from .actions import Action
import numpy
import os
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    lib_file = importer.lazy_import_as('..lib_file', 'lib_file')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')


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
        line = next(Ofile)

        self.numF = int(line)

        while True:
            try:
                line = next(Ofile)
            except StopIteration:
                break

            words = line.split()

            self.state_list.append({'name':words[0],
                'Om':float(words[1]),
                'OmFrag':list(map(float, words[2:]))})
        self.numSt = len(self.state_list)

    def Om_bar_input(self):
        self.read_str("Name of the file with the tden information", "tdenfile", "tden_summ.txt", autocomp=True)

        if self.numSt < 20:
            defw = 7.
        else:
            defw = 15.
        self.read_float("Width of the plot (cm)", 'width', defw)

        self.comps = []
        print("Please enter the different excitation components to be plotted")
        print("    - leave empty to finish")

        colors = ['blue', 'red', 'G', 'O', 'Y']
        for icomp in range(1,1000):
            rstr = self.ret_str('Name of component %i (e.g. MLCT or A-B)'%icomp)
            if rstr == '':
                print(" ... component input finished.")
                break

            color = self.ret_str('Color for plotting', colors[(icomp-1)%5])

            print("\n *** Fragment pairs belonging to %s ***"%rstr)
            print("  Enter two indices between 1 and %i, separated by spaces"%self.numF)
            print("  Leave empty to finish")
            cols = [] # Columns in the OmFrag.txt file
            for jpair in range(1,1000):
                ehstr = self.ret_str("Hole/electron indices for pair %i"%jpair)
                if ehstr == '':
                    print(" ... switching to next component.")
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
            f.write(" %.6f"%(1.-state['Om']))
            f.write("\n")
        f.close()

    def make_tex(self):
        lfile = lib_file.latexfile('Om_bars.tex')
        lfile.write(self.pre())

        # Graphs with energies
        lfile.write("""\\begin{axis}[
      anchor=south west, at={(0,0)},
      height=3.0cm,
      xlabel={}, ylabel={Energy (eV)},
      xtick={-10}, xticklabels={},
      %colorbar % uncomment for colorbar as legend
      ]\n""")

        lfile.write("\\addplot[mark=-, only marks, mark options={scale=2, line width=2pt}]  table[x index=0, y index=2] {Om_bar_data.txt};\n\n")
        lfile.write("% Uncomment to plot oscillator strengths as shading\n")
        lfile.write("%\\addplot+[solid, mesh, point meta=explicit, no markers, line width=5cm, shader=flat corner] table[x expr=\\thisrowno{0}-0.5, y expr=2, meta expr=\\thisrowno{3}] {Om_bar_data.txt};\n")
        lfile.write('\end{axis}\n')

        # Bar graphs with state characters
        lfile.write("""\\begin{axis}[
      anchor=north west, at={(0,0)},
      xlabel={State}, ylabel={Character},
      ytick={0.0,0.2,...,0.8},
      xtick pos=left, xtick align=outside,\n""")
        lfile.write("xtick={1, %i, %i},\n"%(self.numSt//2, self.numSt))
        #lfile.write("xticklables={1, %i, %i},"%(self.numSt//2, self.numSt))
        labs = (self.state_list[0]['name'], self.state_list[self.numSt//2-1]['name'], self.state_list[self.numSt-1]['name'])
        lfile.write("xticklabels={$%s$, $%s$, $%s$},\n"%labs)
        lfile.write("ymin=-0.0, ymax=1,\n")
        lfile.write("ybar stacked, bar width=%.3f cm]\n"%(self['width']/self.numSt/1.5))

        for icomp, comp in enumerate(self.comps):
            lfile.write("\\addplot[ybar, draw=black, fill=%s] table[x index=0, y index=%i] {Om_bar_data.txt};\n"%(comp['color'], icomp+4))
            lfile.write("\\addlegendentry{%s};\n"%comp['name'])
        lfile.write('\n% Uncomment for double excitations\n% ')
        lfile.write("\\addplot[ybar, draw=black, fill=yellow!50] table[x index=0, y index=%i] {Om_bar_data.txt};\n"%(len(self.comps)+4))
        lfile.write('% \\addlegendentry{2-el.};\n')

        lfile.write('''% Handles for arrows
\\node at (axis cs: 1,-0.15) (h1) {};\n''')
        lfile.write("\\node at (axis cs: %i,-0.15) (h2) {};\n"%(self.numSt//2))
        lfile.write("\\node at (axis cs: %i,-0.15) (h3) {};\n"%(self.numSt))
        lfile.write('''\end{axis}\n\n

% Template for drawing orbital pictures
%    You can use the \incMO command defined above
% \\node[orb, anchor=north west, at={(-1cm,-4.8cm)}] (S1e) {incMO};
% \\node[orb, below of=S1e, anchor=north, node distance=2cm] (S1h) {incMO};
% \\node[orb, right of=S1e, anchor=west, xshift=0.7cm, yshift=-0.3cm] (S2e) {incMO};
% \\node[orb, below of=S2e, anchor=north, node distance=2cm] (S2h) {incMO};
% \\node[orb, right of=S2e, anchor=west, xshift=0.7cm, yshift=0.3cm] (S3e) {incMO};
% \\node[orb, below of=S3e, anchor=north, node distance=2cm] (S3h) {incMO};
%

% \draw (S1e) to (S1h);
% \draw (S2e) to (S2h);
% \draw (S3e) to (S3h);

% \draw[arrow] (S1e) to [out=90, in=-90] (h1);
% \draw[arrow] (S2e) to [out=90, in=-90] (h2);
% \draw[arrow] (S3e) to [out=90, in=-90] (h3);
''')

        # Finish
        lfile.write('\end{tikzpicture}\n')
        lfile.post(lvprt=1)
        print("  -> Create plots using: pdflatex %s"%lfile.name)

    def pre(self):
        """
        Code for beginning of tex file.
        """
        str = """\documentclass{standalone}
\\usepackage{xcolor}
\\usepackage{tikz, pgfplots}
\\usetikzlibrary{plotmarks, arrows}
\\usetikzlibrary{calc, decorations, plotmarks, fit, positioning, shapes.geometric}
\pgfplotsset{compat=1.4}

% ===========================================================================
% colour definitions

\definecolor{BL}{HTML}{0099e6}
\definecolor{B}{HTML}{006699}
\definecolor{BD}{HTML}{004466}
\definecolor{BR}{HTML}{0080FF}  % like RoyalBlue

\definecolor{RL}{HTML}{d65454}
\definecolor{R}{HTML}{893636}
\definecolor{RD}{HTML}{ad1737}  % like Maroon

\definecolor{GL}{HTML}{10d05d}
\definecolor{G}{HTML}{009933}
\definecolor{GD}{HTML}{004d1a}

\definecolor{Y}{HTML}{D0D030}
\definecolor{O}{HTML}{FFA500}
\definecolor{T}{HTML}{4DC0FF}
\definecolor{OL}{HTML}{FFAA80}
\definecolor{P}{HTML}{DD88BB}
\definecolor{V}{HTML}{703676}

% ===========================================================================

% settings for the coordinate system
\pgfplotsset{
every linear axis/.append style={
"""
        str += 'width=%.3f cm, height=3.0cm,'%self['width']
        str += "xmin=0, xmax=%i,\n"%(self.numSt+1)
        str += """scale only axis, axis on top,
  yticklabel style={
    /pgf/number format/fixed,
    /pgf/number format/zerofill,
    /pgf/number format/precision=1},
  legend style={at={(1.01,0.00)}, anchor=south west, draw=none, fill=black!5, inner sep=0.5pt, outer sep=0.5pt},
  legend columns=1, legend cell align=left},
}
\pgfplotsset{colormap={CI}{color=(white); color=(blue!50); color=(O!70);}}

% ===========================================================================

\\tikzstyle{orb} = [ellipse, draw, color=black, fill=black!5, inner sep=-1pt, align=center]
\\tikzstyle{arrow} = [-open triangle 45, black, thick]
\\newcommand{\\incMO}[1]{\\includegraphics[width=4.5cm, trim=1cm 1cm 1cm 1cm, clip=true]{#1}}

\\begin{document}
\\begin{tikzpicture}
"""
        return str

class PlotOmBars(Action):
    name = 'plot_om_bars'

    _colt_description = 'Plot Omega matrices as bar graphs'

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
            '..lib_file': 'lib_file',
            '..error_handler': 'error_handler',
    })

    def run():
        theo_header.print_header(__class__._colt_description)
        opt = Om_bar_options('Om_bars.in')
        opt.read_str("Name of the file with the Omega matrix entries", "OmFfile", "OmFrag.txt", autocomp=True)
        opt.read_OmFrag(opt['OmFfile'])
        opt.Om_bar_input()
        opt.flush()
        opt.Om_bar_data()
        opt.make_tex()
