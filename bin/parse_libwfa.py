#!/usr/bin/env python2
"""
Script for parsing libwfa output.
"""

import theo_header, dens_ana_base, input_options
import sys

def ihelp():
    print " parse_libwfa.py <logfile> <type>\n"
    print "  type: qcadc, qctddft, rassi\n"
    print " Command line options:"
    print "  -h, -H, -help: print this help"
    print "  -ifile, -f [dens_ana.in]: name of the input file"
    exit(0)

#--------------------------------------------------------------------------#        
# Input options
#--------------------------------------------------------------------------# 

tmp = sys.argv.pop(0)
#if len(sys.argv) == 0: ihelp()

args2 = []
ifile = 'dens_ana.in'
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "-help"]:
        ihelp()
    elif arg == '-ifile' or arg == '-f':
        ifile = sys.argv.pop(0)
    else:
        args2.append(arg)

ioptions = input_options.libwfa_parse_options(ifile, check_init=False)

if len(args2) >= 1:
    ioptions['rfile'] = args2[0]
    if len(args2) >= 2:
        ioptions['rtype'] = args2[1]

if (not 'rfile' in ioptions) or (not 'rtype' in ioptions):
    ihelp()

theo_header.print_header('Parse libwfa output', ioptions=ioptions)

#--------------------------------------------------------------------------#        
# Parsing and computations
#--------------------------------------------------------------------------#

dena = dens_ana_base.dens_ana_base(ioptions)
#sdena.read_mos()
dena.read_dens()

dena.print_summary()
