#!/usr/bin/env python3
from __future__ import print_function, division
"""
Driver script for analyzing a set of NO files.
"""

import os, sys

from theodore import theo_header, lib_sden, input_options

theo_header.print_header('State density matrix analysis')

def ihelp():
    print(" analyze_sden.py\n")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -ifile, -f  [dens_ana.in]: name of the input file")
    exit(0)

#--------------------------------------------------------------------------#
# Input options
#--------------------------------------------------------------------------#

tmp = sys.argv.pop(0)

ifile = 'dens_ana.in'
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "--help"]:
        ihelp()
    elif arg == '-ifile' or arg == '-f':
        ifile = sys.argv.pop(0)

if not os.path.exists(ifile):
    print('Input file %s not found!'%ifile)
    print('Please create this file using theoinp or specify its location using -ifile\n')
    ihelp()

ioptions = input_options.sden_ana_options(ifile)

#--------------------------------------------------------------------------#
# Parsing and computations
#--------------------------------------------------------------------------#

sdena = lib_sden.sden_ana(ioptions)
if 'mo_file' in ioptions: sdena.read_mos()
sdena.read_dens()

if ioptions['NO_ana']:  sdena.compute_all_NO()
if ioptions['AD_ana']:  sdena.compute_all_AD()
if ioptions['BO_ana']:  sdena.compute_all_BO()

if ioptions['comp_rho']: sdena.compute_rho()

#--------------------------------------------------------------------------#
# Print out
#--------------------------------------------------------------------------#
if ioptions['pop_ana']: sdena.print_all_pop_table()
if ioptions['BO_ana']:  sdena.print_all_BO()
if ioptions['mo_pop_type'] > 0: sdena.print_mo_pops(ioptions['mo_pop_type'])

sdena.print_summary()
