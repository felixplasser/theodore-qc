#!/usr/bin/env python3
"""
Script for parsing libwfa output.
"""
from __future__ import print_function, division
import sys

from .. import theo_header, dens_ana_base, input_options, error_handler
from .actions import Action

def ihelp():
    print(" parse_libwfa.py <logfile> <type>\n")
    print("  type: qcadc, qctddft, qctda, rassi\n")
    print(" Command line options:")
    print("  -h, -H, -help: print this help")
    print("  -ifile, -f [dens_ana.in]: name of the input file")
    exit(0)

class ParseLibwfa(Action):    

    name = 'parse_libwfa'

    _questions = """
    logfile = :: existing_file
    typ    = :: str, optional :: qcadc, qctddft, qctda, rassi
    ifile = dens_ana.in :: existing_file, alias=f
    """

    def run(logfile, typ, ifile):
        #--------------------------------------------------------------------------#
        # Input options
        #--------------------------------------------------------------------------#
        
        ioptions = input_options.libwfa_parse_options(ifile, check_init=False)
        
        ioptions['rfile'] = logfile
        if typ is not None:
            ioptions['rtype'] = typ
        
        if ioptions['rtype'] == 'qctda':
            ioptions['TDA'] = True
            ioptions['rtype'] = 'qctddft'
        
        if (not 'rfile' in ioptions) or (not 'rtype' in ioptions):
            ihelp()
        
        theo_header.print_header('Parse libwfa output', ioptions=ioptions)
        
        #--------------------------------------------------------------------------#
        # Parsing and computations
        #--------------------------------------------------------------------------#
        
        dena = dens_ana_base.dens_ana_base(ioptions)
        #sdena.read_mos()
        
        dena.read_dens()
        
        dena.print_summary()
