#!/usr/bin/env python3

"""
Unit tests for TheoDORE.
Activate this by running pytest [-v/-s] in the EXAMPLES directory.
"""

import pytest
import os

def test_header():
    from theodore import theo_header
    
    # TODO: make sure this is using the correct version here
    print
    #theo_header.print_header('Unit tests')
    print("THEODIR: ", os.environ["THEODIR"])

class Test1:
    def test_import(self):
        """
        Test import statements.
        """
        from theodore import atominfo
        from theodore import cclib_interface
        from theodore import dens_ana_base
        from theodore import error_handler
        from theodore import fchk_parser
        from theodore import input_options
        from theodore import lib_diab
        from theodore import lib_exciton
        from theodore import lib_file
        from theodore import lib_mo
        from theodore import lib_NICS
        from theodore import lib_plot
        from theodore import lib_sden
        from theodore import lib_soc
        from theodore import lib_struc
        from theodore import lib_tden
        from theodore import lib_util
        from theodore import OB_repl
        from theodore import Om_descriptors
        from theodore import orbkit_interface
        from theodore import pop_ana
        from theodore import theo_header
        from theodore import units
