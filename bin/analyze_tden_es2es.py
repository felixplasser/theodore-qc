#!/usr/bin/env python3
"""
Transition density matrix analysis between two excited states.

This is computed analogously to the electron/hole densities:
D^IJ = (D^0I)^T * D^0J - D^0I * (D^0J)^T
"""

import os, sys

from theodore import theo_header, lib_tden, lib_exciton, input_options, error_handler

def ihelp():
    print(" analyze_tden_es2es.py")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -ifile, -f [dens_ana.in]: name of the input file")
    print("  -r [1]: index of reference state for state-to-state analysis")
    exit(0)

iref = 1
ifile = 'dens_ana.in'

arg=sys.argv.pop(0)
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "--help"]:
        ihelp()
    elif arg == '-ifile' or arg == '-f':
        ifile = sys.argv.pop(0)
    elif arg == '-r':
        iref = int(sys.argv.pop(0))
    else:
        raise error_handler.ElseError(arg, 'command line option')

if not os.path.exists(ifile):
    print('Input file %s not found!'%ifile)
    print('Please create this file using theoinp or specify its location using -ifile\n')
    ihelp()

ioptions = input_options.tden_ana_options(ifile)
theo_header.print_header('Trans. dens. mat. analysis (excited-to-excited)', ioptions=ioptions, cfile=__file__)

tdena = lib_tden.tden_ana(ioptions)
if 'mo_file' in ioptions: tdena.read_mos()

tdena.read_dens()
tdena.compute_es2es_tden(iref=iref)

if 'at_lists' in ioptions or ioptions['eh_pop'] >= 1:
    tdena.compute_all_OmAt()

if 'at_lists' in ioptions:
    tdena.compute_all_OmFrag()
    if ioptions['print_OmFrag']: tdena.fprint_OmFrag()

if ioptions['comp_ntos']:  tdena.compute_all_NTO()
if ioptions['comp_p_h_dens']: tdena.compute_p_h_dens()
if ioptions['comp_rho0n']: tdena.compute_rho_0_n()
if 'Phe' in ioptions['prop_list']:
    tdena.compute_all_Phe()

#--------------------------------------------------------------------------#
# Print-out
#--------------------------------------------------------------------------#    
tdena.print_all_eh_pop()

tdena.print_summary()