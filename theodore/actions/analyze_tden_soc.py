"""
Driver script for transition density matrix analysis.
"""
from __future__ import print_function, division
import os, sys, time

from .actions import Action
from .theotools import timeit

from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_soc = importer.lazy_import_as('..lib_soc', 'lib_soc')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    input_options = importer.lazy_import_as('..input_options', 'input_options')


#--------------------------------------------------------------------------#
# Parsing and computations
#--------------------------------------------------------------------------#

class AnalyzeTdenSoc(Action):
    name = 'analyze_tden_soc'

    _colt_description = '1TDM analysis for spin-orbit coupled states'

    _user_input = """
    ifile = dens_ana.in :: existing_file, alias=f
    spin_comp = False :: bool, alias=s
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_soc': 'lib_soc',
            '..error_handler': 'error_handler',
            '..input_options': 'input_options',
    })

    @timeit
    def run(ifile, spin_comp):
        ioptions = input_options.tden_ana_options(ifile)
        theo_header.print_header(__class__._colt_description, ioptions=ioptions, cfile=__file__)
        
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
