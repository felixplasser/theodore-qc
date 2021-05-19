#!/usr/bin/env python3
"""
Script for parsing libwfa output.
"""
from __future__ import print_function, division
import sys

from .. import theo_header, dens_ana_base, input_options, error_handler

def ihelp():
    print(" parse_libwfa.py <logfile> <type>\n")
    print("  type: qcadc, qctddft, qctda, rassi\n")
    print(" Command line options:")
    print("  -h, -H, -help: print this help")
    print("  -ifile, -f [dens_ana.in]: name of the input file")
    exit(0)

    
def parse_libwfa():
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
    
    if ioptions['rtype'] == 'qctda':
        ioptions['TDA'] = True
        ioptions['rtype'] = 'qctddft'
    
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
