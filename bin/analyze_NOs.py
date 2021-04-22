"""
Driver script for analyzing a set of NO files.
"""

from __future__ import print_function, division
import sys

from .. import theo_header, lib_sden, input_options

import argparse


def get_commandline_args():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('mo_file', nargs='+', help='')
    parser.add_argument('-ifile', default='dens_ana.in', help='name of the input file')
    args = parser.parse_args()
    return args.mo_file, args.ifile


def analyze_nos():
    theo_header.print_header('NO file analysis')
    no_files, ifile = get_commandline_args()
    # set options
    ioptions = input_options.sden_ana_options(ifile, check_init=False)
    ioptions['rtype'] = 'nos'
    ioptions['mo_file'] = no_files[0]
    ioptions['ana_files'] = no_files[1:]

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
