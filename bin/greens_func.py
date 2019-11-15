#!/usr/bin/env python3
from __future__ import print_function, division
"""
Driver script for analyzing a set of NO files.
"""

import os, sys

from theodore import theo_header, lib_green

theo_header.print_header('Construction of Green\'s functions')

def ihelp():
    print(" greens_func.py\n")
    print(" Command line options:")
    print("  -h, -H, -help: print this help")
    print("  -ifile, -f  [dens_ana.in]: name of the input file")
    print("  -E: energy for Green's function")
    exit(0)

#--------------------------------------------------------------------------#
# Input options
#--------------------------------------------------------------------------#

tmp = sys.argv.pop(0)
EG = None

ifile = 'dens_ana.in'
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "-help"]:
        ihelp()
    elif arg == '-ifile' or arg == '-f':
        ifile = sys.argv.pop(0)
    elif arg == '-E':
        EG = float(sys.argv.pop(0))

if not os.path.exists(ifile):
    print('Input file %s not found!'%ifile)
    print('Please create this file using theoinp or specify its location using -ifile\n')
    ihelp()

ioptions = lib_green.green_options(ifile)
ioptions['rtype'] = 'nos'
ioptions['Om_formula'] = 2

#--------------------------------------------------------------------------#
# Parsing and computations
#--------------------------------------------------------------------------#

ga = lib_green.green_ana(ioptions)
ga.read_mos()


energies = [0.5*(-0.252)+0.5*0.049, 0.75*(-0.252)+0.25*0.049]

ga.compute_G(energies)
