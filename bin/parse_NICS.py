#!/usr/bin/env python3
"""
Driver script for parsing NICS values and various plotting routines.
"""

from theodore import theo_header, lib_NICS

logfile = 'gaussian.log'

nv = lib_NICS.NICS_parser_g09(logfile)
nv.print_data()