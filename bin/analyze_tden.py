#!/usr/bin/env python3
"""
Driver script for transition density matrix analysis.
"""
from __future__ import print_function, division
import os, sys, time

# Python 2/3 compatibility
try:
    from time import process_time
except ImportError:
    from time import clock as process_time
try:
    from time import perf_counter
except ImportError:
    from time import time as perf_counter

from theodore import theo_header, lib_tden, lib_exciton, input_options, error_handler

(tc, tt) = (process_time(), perf_counter())

def ihelp():
    print(" analyze_tden.py")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -ifile, -f [dens_ana.in]: name of the input file")
    exit(0)

#--------------------------------------------------------------------------#
# Parsing and computations
#--------------------------------------------------------------------------#

ifile = 'dens_ana.in'

arg=sys.argv.pop(0)
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "--help"]:
        ihelp()
    elif arg == '-ifile' or arg == '-f':
        ifile = sys.argv.pop(0)
    else:
        raise error_handler.ElseError(arg, 'command line option')

if not os.path.exists(ifile):
    print('Input file %s not found!'%ifile)
    print('Please create this file using theoinp or specify its location using -ifile\n')
    ihelp()

ioptions = input_options.tden_ana_options(ifile)
theo_header.print_header('Transition density matrix analysis', ioptions=ioptions)

tdena = lib_tden.tden_ana(ioptions)
if 'mo_file' in ioptions: tdena.read_mos()

tdena.read_dens()

if 'at_lists' in ioptions or ioptions['eh_pop'] >= 1:
    tdena.compute_all_OmAt()

if 'at_lists' in ioptions:
    tdena.compute_all_OmFrag()
    if ioptions['print_OmFrag']: tdena.fprint_OmFrag()

if ioptions['comp_ntos']:  tdena.compute_all_NTO()
if ioptions['comp_dntos']: tdena.compute_all_DNTO()
if ioptions['comp_p_h_dens']: tdena.compute_p_h_dens()
if ioptions['comp_rho0n']: tdena.compute_rho_0_n()

if 'RMSeh' in ioptions.get('prop_list') or 'MAeh' in ioptions.get('prop_list') or 'Eb' in ioptions.get('prop_list'):
    exca = lib_exciton.exciton_analysis()
    exca.get_distance_matrix(tdena.struc)
    tdena.analyze_excitons(exca)

if 'Phe' in ioptions['prop_list']:
    tdena.compute_all_Phe()

#--------------------------------------------------------------------------#
# Print-out
#--------------------------------------------------------------------------#

tdena.print_all_eh_pop()

tdena.print_summary()

#print 'Finished at ' + time.asctime()

print("CPU time: % .1f s, wall time: %.1f s"%(process_time() - tc, perf_counter() - tt))
