"""
Convert a log-file to Molden format with the help of cclib.
"""
from __future__ import print_function, division
import sys

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    cclib_interface = importer.lazy_import_as('..cclib_interface', 'cclib_interface')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    lib_mo = importer.lazy_import_as('..lib_mo', 'lib_mo')
    input_options = importer.lazy_import_as('..input_options', 'input_options')


class CC2Molden(Action):

    name = 'cc2molden'

    _user_input = """
    logfile = :: existing_file
    """
    _colt_description = "Convert a log-file to Molden format with the help of cclib."

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..cclib_interface': 'cclib_interface',
            '..error_handler': 'error_handler',
            '..input_options': 'input_options',
            '..lib_mo': 'lib_mo',
    })

    def run(logfile):
        theo_header.print_header(__class__._colt_description)
        print("  WARNING: This script is not well-tested and might fail for some of the quantum chemistry codes.")
        
        ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
        ioptions['rtype'] = 'cclib'
        ioptions['rfile'] = logfile
        
        ccparser = cclib_interface.file_parser_cclib(ioptions)
        
        errcode = ccparser.check()
        if errcode == 0:
            mos = ccparser.read_mos()
            mos.write_molden_file(fname="cc.mld")
            print("\n Finished: molden format file cc.mld written.")
        else:
            print(" Conversion to Molden format not possible!")
            print(" %s does not contain all required information"%logfile)
