#!/usr/bin/env python3
"""
Script for parsing libwfa output.
"""
from __future__ import print_function, division
import sys

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    dens_ana_base = importer.lazy_import_as('..dens_ana_base', 'dens_ana_base')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    input_options = importer.lazy_import_as('..input_options', 'input_options')

class ParseLibwfa(Action):    

    name = 'parse_libwfa'

    _colt_description = 'Parse libwfa output from Q-Chem or OpenMolcas'

    _user_input = """
    # Logfile from Q-Chem or OpenMolcas
    logfile = :: existing_file
    # Type of calculation (qcadc, qctddft, qctda, rassi)
    typ    = :: str :: qcadc, qctddft, qctda, rassi
    # Input file
    ifile = :: file, optional, alias=f
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..dens_ana_base': 'dens_ana_base',
            '..error_handler': 'error_handler',
            '..input_options': 'input_options',
    })

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

        theo_header.print_header(__class__._colt_description, ioptions=ioptions)
        
        #--------------------------------------------------------------------------#
        # Parsing and computations
        #--------------------------------------------------------------------------#
        
        dena = dens_ana_base.dens_ana_base(ioptions)
        #sdena.read_mos()
        
        dena.read_dens()
        
        dena.print_summary()
