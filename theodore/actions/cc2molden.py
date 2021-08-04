"""
Convert a log-file to Molden format with the help of cclib.
"""
from __future__ import print_function, division
import sys
import argparse

from .. import theo_header, cclib_interface, input_options, error_handler, lib_mo
from .actions import Action

class CC2Molden(Action):

    name = 'cc2molden'
    _questions = """
    logfile = :: existing_file
    """
    _colt_description = "Convert a log-file to Molden format with the help of cclib."

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
