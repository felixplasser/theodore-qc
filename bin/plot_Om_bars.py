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
    ## \param fname file with the data produced in a previous analyze_tden.py run
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
            print "  Enter indices between 1 and %i, leave empty to finish"%self.numF
            cols = [] # Columns in the OmFrag.txt file
            for jpair in range(1,1000):
                jh = self.ret_str("Hole index for pair %i"%jpair)
                if jh == '':
                    print " ... switching to next pair."
                    break
                je = self.ret_str("Electron index for pair %i"%jpair)
                cols.append((int(jh)-1) * self.numF + (int(je)-1))

            self.comps.append({'name':rstr, 'color':color, 'cols':cols})
        print self.comps

    def Om_bar_data(self):
        """
        Create file with the data to be plotted.
        """
        f = open('Om_bar_data.txt', 'w')
        ist = 1
        for state in self.state_list:
            f.write("%3i"%ist)
            ist += 1
            for comp in self.comps:
                wt = sum(state['OmFrag'][icol] for icol in comp['cols'])
                f.write(" %.6f"%wt)
            f.write("\n")
        f.close()

    def make_tex(self):
        lfile = lib_file.latexfile('Om_bars.tex')
        lfile.write(self.pre())
        lfile.write("""\\begin{axis}[
      anchor=north west, at={(0,0)},
      xlabel={State}, ylabel={Character},
      xticklabels={}, ytick={0.0,0.2,...,0.8},
""")
        lfile.write("xmin=0, xmax=%i,\n"%(self.numSt+1))
        lfile.write("ymin=-0.0, ymax=1, xtick pos=both,\n")
        lfile.write("ybar stacked, bar width=%.3f cm]\n"%(self['width']/self.numSt/1.2))

        for icomp, comp in enumerate(self.comps):
            lfile.write("\\addplot[ybar, draw=none, fill=%s] table[x index=0, y index=%i] {Om_bar_data.txt};\n"%(comp['color'], icomp+1))
            lfile.write("\\addlegendentry{%s};\n"%comp['name'])
        lfile.write('\end{axis}\n')

        lfile.write('\end{tikzpicture}\n')
        lfile.post(lvprt=1)
        print "  -> Create plots using: pdflatex %s"%lfile.name

    def pre(self):
        """
        Code for beginning of tex file.
        """
        str = """\documentclass[10pt]{article}
\usepackage{xcolor}

% here comes the graphics package with a number of useful libraries
\usepackage{tikz, pgfplots}
\usetikzlibrary{calc, decorations, plotmarks, fit, positioning, shapes.geometric}
\pgfplotsset{compat=1.4}

% setup preview package
\usepackage[active,tightpage]{preview}
\PreviewEnvironment{tikzpicture}                % extract tikzpicture environments
\setlength\PreviewBorder{0pt}                   % make tightly fitting pages around

% ===================================================================================
% ===================================================================================

% settings for the coordinate system
\pgfplotsset{
every linear axis/.append style={
"""
        str += 'width=%.3f cm, height=2.0cm,'%self['width']
        str += """scale only axis,                  % given size applies to axis, not to whole plot
  axis on top,
%
  yticklabel style={                %
    /pgf/number format/fixed,       %
    /pgf/number format/zerofill,    %
    /pgf/number format/precision=1  %
  },
  legend style={at={(1.01,0.00)}, anchor=south west, draw=none, fill=black!5, inner sep=0.5pt, outer sep=0.5pt},
  legend columns=1,
  legend cell align=left,                   % left-align legend entries
  solid,                                            % with lines
  clip mode=individual,
  line join=bevel,
},
every mark/.append style={very thin}
}

\\begin{document}
\\begin{tikzpicture}
"""
        return str

def run_plot():
    opt = Om_bar_options('plot.in')
    opt.read_OmFrag()
    opt.Om_bar_input()
    opt.Om_bar_data()
    opt.make_tex()


if __name__ == '__main__':
    theo_header.print_header('Plot Omega matrices as bar graphs')

    run_plot()
