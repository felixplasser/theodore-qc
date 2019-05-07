#!/usr/bin/env python3
"""
Driver script for transition density matrix analysis.
"""
from __future__ import print_function, division
import os, sys, time

from theodore import theo_header, lib_soc, input_options, error_handler

(tc, tt) = (time.clock(), time.time())

def ihelp():
    print(" analyze_tden_soc.py")
    print(" Command line options:")
    print("  -h, -H, -help: print this help")
    print("  -ifile, -f [dens_ana.in]: name of the input file")
    print("  -s: analyze spin components")
    exit(0)

#--------------------------------------------------------------------------#
# Parsing and computations
#--------------------------------------------------------------------------#

ifile = 'dens_ana.in'
spin_comp = False

arg=sys.argv.pop(0)
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "-help"]:
        ihelp()
    elif arg == '-ifile' or arg == '-f':
        ifile = sys.argv.pop(0)
    elif arg == '-s':
        spin_comp = True
    else:
        raise error_handler.ElseError(arg, 'command line option')

if not os.path.exists(ifile):
    print('Input file %s not found!'%ifile)
    print('Please create this file using theoinp or specify its location using -ifile\n')
    ihelp()

ioptions = input_options.tden_ana_options(ifile)
theo_header.print_header('1TDM analysis for spin-orbit coupled states', ioptions=ioptions)

tdena = lib_soc.tden_ana_soc(ioptions)
if 'mo_file' in ioptions: tdena.read_mos()
tdena.read_dens()
tdena.compute_all_OmAt(fullmat=True)
tdena.soc_transform()

tdena.print_info('mch')
if spin_comp:
    tdena.print_info('aa')
    tdena.print_info('bb')
    tdena.print_info('ab')
    tdena.print_info('ba')
tdena.print_info('soc')

#print 'Finished at ' + time.asctime()

print("CPU time: % .1f s, wall time: %.1f s"%(time.clock() - tc, time.time() - tt))
