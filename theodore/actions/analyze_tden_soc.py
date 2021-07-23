"""
Driver script for transition density matrix analysis.
"""
from __future__ import print_function, division
import os, sys, time

from .. import theo_header, lib_soc, input_options, error_handler
from .actions import Action
from .theotools import timeit



def get_commandline_args():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('-ifile', default='dens_ana.in', help='name of the input file')
    parser.add_argument('-s', action='store_true', help='')
    args = parser.parse_args()
    return isfile(args.ifile), args.s


#--------------------------------------------------------------------------#
# Parsing and computations
#--------------------------------------------------------------------------#

class AnalyzeTdenSoc(Action):
    name = 'analyze_tden_soc'

    _colt_description = '1TDM analysis for spin-orbit coupled states'

    _questions = """
    ifile = dens_ana.in :: existing_file, alias=f
    spin_comp = False :: bool, alias=s
    """

    @timeit
    def run(ifile, spin_comp):
        ioptions = input_options.tden_ana_options(ifile)
        theo_header.print_header(_colt_description, ioptions=ioptions, cfile=__file__)
        
        tdena = lib_soc.tden_ana_soc(ioptions)
        if 'mo_file' in ioptions: tdena.read_mos()
        tdena.read_dens()
        tdena.compute_all_OmAt(fullmat=True)
        tdena.soc_transform()
        
        tdena.print_info('mch')
        if spin_comp is True:
            tdena.print_info('aa')
            tdena.print_info('bb')
            tdena.print_info('ab')
            tdena.print_info('ba')
        tdena.print_info('soc')
