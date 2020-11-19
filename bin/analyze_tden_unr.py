#!/usr/bin/env python3
"""
Driver script for transition density matrix analysis in the case of unrestricted orbitals.
Authors: Felix Plasser, Sebastian Mai
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
    print(" analyze_tden_unr.py")
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
theo_header.print_header('Transition density matrix analysis (UHF/UKS)', ioptions=ioptions)

ioptions['jmol_orbitals'] = False

# ALPHA spin
print("\nRunning alpha-spin analysis in directory ALPHA")
ioptions['spin'] = 1

tdena_alpha = lib_tden.tden_ana(ioptions)
if 'mo_file' in ioptions: tdena_alpha.read_mos(spin=1)
tdena_alpha.read_dens()

try:
    os.mkdir('ALPHA')
except OSError:
    pass
os.chdir('ALPHA')
if 'at_lists' in ioptions:
    tdena_alpha.compute_all_OmFrag()
    if ioptions['print_OmFrag']: tdena_alpha.fprint_OmFrag()
if ioptions['comp_ntos']:  tdena_alpha.compute_all_NTO()

print("\n *** ALPHA-spin results ***")
tdena_alpha.print_summary()
os.chdir('..')

# BETA spin
print("\nRunning beta-spin analysis in directory BETA")
ioptions['spin'] = -1

tdena_beta = lib_tden.tden_ana(ioptions)
if 'mo_file' in ioptions: tdena_beta.read_mos(spin=-1)
tdena_beta.read_dens()

try:
    os.mkdir('BETA')
except OSError:
    pass
os.chdir('BETA')
if 'at_lists' in ioptions:
    tdena_beta.compute_all_OmFrag()
    if ioptions['print_OmFrag']: tdena_beta.fprint_OmFrag()
if ioptions['comp_ntos']:  tdena_beta.compute_all_NTO()

print("\n *** BETA-spin results ***")
tdena_beta.print_summary()
os.chdir('..')

# ALPHA+BETA
print("Starting spin-summed analysis")
# Add the alpha values on top of the beta values
for i, state in enumerate(tdena_beta.state_list):
    # Add the things that are additive
    for aprop in ['Om', 'OmAt', 'OmFrag', 'S_HE']:
        if aprop in state:
            state[aprop] += tdena_alpha.state_list[i][aprop]
    try:
        state['Z_HE'] = 2.**(state['S_HE'])
    except KeyError:
        pass

    # Delete the things that are non-additive
    for dprop in ['tden', 'PRNTO', 'Om_desc']:
        if dprop in state:
            del state[dprop]

if 'at_lists' in ioptions:
    tdena_beta.compute_all_OmFrag()
    if ioptions['print_OmFrag']: tdena_beta.fprint_OmFrag()

print("\n *** Spin-summed results ***")
tdena_beta.print_summary()

print("CPU time: % .1f s, wall time: %.1f s"%(process_time() - tc, perf_counter() - tt))
