#!/usr/bin/env python3
from __future__ import print_function, division
"""
Driver script for analyzing a set of NO files.
"""

import os, sys

from theodore import theo_header, lib_green, units

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

ihomo = ga.mos.ret_ihomo()
ehomo = ga.mos.ens[ihomo]
elumo = ga.mos.ens[ihomo+1]

# Fermi energy
EF = 0.5 * ehomo + 0.5 * elumo

print(ehomo, elumo, EF)

#energies = [(1.-x)*ehomo + x*elumo for x in [0.1, 0.25, 0.5, 0.75, 0.9]]

energies = [EF + x / units.energy['eV'] for x in [-1, -0.5, -0.25, 0, 0.25, 0.5, 1]]

ga.compute_G(energies)
