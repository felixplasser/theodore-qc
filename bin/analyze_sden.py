#!/usr/bin/python
"""
Driver script for analyzing a set of NO files.
"""

import theo_header, lib_sden, input_options
import os, sys

theo_header.print_header('NO file analysis')

def ihelp():
    print " analyze_sden.py\n"
    print " Command line options:"
    print "  -h, -H, -help: print this help"
    print "  -ifile [dens_ana.in]: name of the input file"
    exit(0)

#--------------------------------------------------------------------------#        
# Input options
#--------------------------------------------------------------------------# 

tmp = sys.argv.pop(0)

ifile = 'dens_ana.in'
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "-help"]:
        ihelp()
    elif arg == '-ifile':
        ifile = sys.argv.pop(0)

if not os.path.exists(ifile):
    print 'Input file %s not found!'%ifile
    print 'Please create this file using theoinp or specify its location using -ifile\n'
    ihelp()
    
ioptions = input_options.sden_ana_options(ifile)

#--------------------------------------------------------------------------#        
# Parsing and computations
#--------------------------------------------------------------------------#

sdena = lib_sden.sden_ana(ioptions)
if 'mo_file' in ioptions: sdena.read_mos()
sdena.read_dens()

#sdena.print_all_mullpop()

if ioptions['AD_ana']:  sdena.compute_all_AD()

#--------------------------------------------------------------------------#        
# Print out
#--------------------------------------------------------------------------#
if ioptions['pop_ana']: sdena.print_all_pop_table()

sdena.print_summary()
