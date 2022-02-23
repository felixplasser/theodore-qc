"""
Driver script for analyzing a set of NO files.
"""

from __future__ import print_function, division
import sys

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_sden = importer.lazy_import_as('..lib_sden', 'lib_sden')
    input_options = importer.lazy_import_as('..input_options', 'input_options')

class AnalyzeNOs(Action):

    name = 'analyze_nos'

    _colt_description = 'Analysis of natural orbital (NO) files'

    _user_input = """
    # List of NO files in Molden format
    no_files = :: list(existing_file)
    # Input file (optional)
    ifile = dens_ana.in :: file, alias=f
    # Reference MO file for computing AO overlap matrix
    ref = :: existing_file, optional, alias=r
    # Multiply occupations with this factor
    occ_fac = :: float, optional, alias=o
    # Use if unrestricted orbitals are present
    unrestricted = false :: bool, alias=u
    # Interpret energies as occupations
    rd_ene = false :: bool, alias=e
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_sden': 'lib_sden',
            '..input_options': 'input_options'
    })

    def run(no_files, ifile, ref, occ_fac, unrestricted, rd_ene):
        theo_header.print_header(__class__._colt_description, cfile=__file__)

        # set options
        ioptions = input_options.sden_ana_options(ifile, check_init=False)
        ioptions['rtype'] = 'nos'
        if not occ_fac is None:
            ioptions['occ_fac'] = occ_fac
        ioptions['unrestricted'] = unrestricted
        ioptions['rd_ene'] = rd_ene

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
        if unrestricted:
            sdena.read_mos(spin=1)
        else:
            sdena.read_mos()
        sdena.read_dens()

        if unrestricted:
            sdena.compute_all_NO()
        if ioptions['AD_ana']:
            sdena.compute_all_AD()
        if ioptions['pop_ana']:
            sdena.print_all_pop_table()
        if ioptions['BO_ana']:
            sdena.compute_all_BO()
            sdena.print_all_BO()

        sdena.print_summary()
