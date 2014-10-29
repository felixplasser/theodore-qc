#!/usr/bin/python
"""
Driver script for analyzing a set of NO files.
"""

import theo_header, lib_sden, input_options
import sys

theo_header.print_header('NO file analysis')

def ihelp():
    print " analyze_NOs.py <NO_file_ref> [<NO_file2> ...]\n"
    print " Command line options:"
    print "  -h, -H, -help: print this help"
    print "  -ifile [dens_ana.in]: name of the input file"
    exit(0)

#--------------------------------------------------------------------------#        
# Input options
#--------------------------------------------------------------------------# 

tmp = sys.argv.pop(0)

no_files = []
ifile = 'dens_ana.in'
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "-help"]:
        ihelp()
    elif arg == '-ifile':
        ifile = sys.argv.pop(0)
    else:
        no_files.append(arg)

if len(no_files) == 0: ihelp()

ioptions = input_options.sden_ana_options(ifile)
ioptions['rtype'] = 'nos'
ioptions['mo_file'] = no_files[0]
ioptions['no_files'] = no_files

#--------------------------------------------------------------------------#        
# Parsing and computations
#--------------------------------------------------------------------------#

sdena = lib_sden.sden_ana(ioptions)
sdena.read_mos()
sdena.read_dens()

#sdena.print_all_mullpop()

if ioptions['pop_ana']: sdena.print_all_pop_table()
if ioptions['AD_ana']:  sdena.compute_all_AD()

sdena.print_summary()
