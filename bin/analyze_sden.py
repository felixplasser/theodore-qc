"""
Driver script for analyzing a set of NO files.
"""

from __future__ import print_function, division
import os, sys
from .. import theo_header, lib_sden, input_options
from .theotools import get_ifile_commandline
import argparse


def analyze_sden():
    # header
    theo_header.print_header('State density matrix analysis')
    ifile = get_ifile_commandline()
    #
    ioptions = input_options.sden_ana_options(ifile)

    #--------------------------------------------------------------------------#
    # Parsing and computations
    #--------------------------------------------------------------------------#

    sdena = lib_sden.sden_ana(ioptions)
    if 'mo_file' in ioptions: sdena.read_mos()
    sdena.read_dens()

    if ioptions['NO_ana']:  sdena.compute_all_NO()
    if ioptions['AD_ana']:  sdena.compute_all_AD()
    if ioptions['BO_ana']:  sdena.compute_all_BO()

    if ioptions['comp_rho']: sdena.compute_rho()

    #--------------------------------------------------------------------------#
    # Print out
    #--------------------------------------------------------------------------#
    if ioptions['pop_ana']: sdena.print_all_pop_table()
    if ioptions['BO_ana']:  sdena.print_all_BO()
    if ioptions['mo_pop_type'] > 0: sdena.print_mo_pops(ioptions['mo_pop_type'])

    sdena.print_summary()
