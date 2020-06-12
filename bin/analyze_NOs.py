#!/usr/bin/env python3
"""
Driver script for analyzing a set of NO files.
"""

from __future__ import print_function, division
import sys

from theodore import theo_header, lib_sden, input_options


theo_header.print_header('NO file analysis')

def ihelp():
    print(" analyze_NOs.py <MO_file> [<NO_file_ref> <NO_file2> ...]\n")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -ifile [dens_ana.in]: name of the input file")
    exit(0)

#--------------------------------------------------------------------------#        
# Input options
#--------------------------------------------------------------------------# 

tmp = sys.argv.pop(0)
if len(sys.argv) == 0: ihelp()

no_files = []
ifile = 'dens_ana.in'
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "--help"]:
        ihelp()
    elif arg == '-ifile':
        ifile = sys.argv.pop(0)
    else:
        no_files.append(arg)

#if len(no_files) == 0: ihelp()

ioptions = input_options.sden_ana_options(ifile, check_init=False)
ioptions['rtype'] = 'nos'
ioptions['mo_file'] = no_files[0]
ioptions['ana_files'] = no_files[1:]

#--------------------------------------------------------------------------#        
# Parsing and computations
#--------------------------------------------------------------------------#

sdena = lib_sden.sden_ana(ioptions)
sdena.read_mos()
sdena.read_dens()

if ioptions['AD_ana']:  sdena.compute_all_AD()
if ioptions['pop_ana']: sdena.print_all_pop_table()
if ioptions['BO_ana']:
    sdena.compute_all_BO()
    sdena.print_all_BO()

sdena.print_summary()
