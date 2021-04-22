#!/usr/bin/env python3
"""
Perform a fragment charge difference analysis, following
  A. A. Voityuk, N. Roesch J. Chem. Phys. 2002, 117, 5607.
"""

from __future__ import print_function, division
from .. import theo_header, input_options,  error_handler, lib_diab
import sys, os, numpy

def ihelp():
    print(" analyze_tden.py")
    print(" Command line options:")
    print("  -h, -H, -help: print this help")
    print("  -ifile, -f [dens_ana.in]: name of the input file")
    exit(0)

def fcd():
    theo_header.print_header('Fragment charge difference analysis')

    ifile = 'fcd.in'

    arg=sys.argv.pop(0)
    while len(sys.argv)>0:
        arg = sys.argv.pop(0)
        if arg in ["-h", "-H", "-help"]:
            ihelp()
        elif arg == '-ifile' or arg == '-f':
            ifile = sys.argv.pop(0)
        else:
            raise error_handler.ElseError(arg, 'command line option')

    if not os.path.exists(ifile):
        print('Input file %s not found!'%ifile)
        print('Please create this file using theoinp or specify its location using -ifile\n')
        ihelp()

    ioptions = input_options.fcd_ana_options(ifile)

    fcda = lib_diab.fcd_ana(ioptions)
    fcda.read_mos()
    fcda.read_dens()

    fcda.do_pop_ana()
    fcda.do_fcd()
