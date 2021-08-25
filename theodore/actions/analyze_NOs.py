"""
Driver script for analyzing a set of NO files.
"""

from __future__ import print_function, division
import sys

from .actions import Action
from .. import theo_header, lib_sden, input_options

class AnalyzeNOs(Action):

    name = 'analyze_nos'

    _colt_description = 'Analysis of natural orbital (NO) files'

    _questions = """
    # List of NO files in Molden format
    no_files = :: list(existing_file)
    # Input file (optional)
    ifile = dens_ana.in :: file, alias=f
    # Reference MO file for computing overlap matrix
    ref = :: existing_file, optional, alias=r
    """

    def run(no_files, ifile, ref):
        theo_header.print_header(__class__._colt_description, cfile=__file__)

        # set options
        ioptions = input_options.sden_ana_options(ifile, check_init=False)
        ioptions['rtype'] = 'nos'

        # optionally use a manually specified MO file for computing the AO overlap matrix
        if ref is None:
            ioptions['mo_file'] = no_files[0]
        else:
            ioptions['mo_file'] = ref
        ioptions['ana_files'] = no_files

        #--------------------------------------------------------------------------#        
        # Parsing and computations
        #--------------------------------------------------------------------------#

        sdena = lib_sden.sden_ana(ioptions)
        sdena.read_mos()
        sdena.read_dens()

        if ioptions['AD_ana']:
            sdena.compute_all_AD()
        if ioptions['pop_ana']:
            sdena.print_all_pop_table()
        if ioptions['BO_ana']:
            sdena.compute_all_BO()
            sdena.print_all_BO()

        sdena.print_summary()
