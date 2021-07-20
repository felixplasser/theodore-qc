#!/usr/bin/env python3
"""
Perform a fragment charge difference analysis, following
  A. A. Voityuk, N. Roesch J. Chem. Phys. 2002, 117, 5607.
"""

from __future__ import print_function, division

from .. import theo_header, input_options,  lib_diab
from .actions import Action


class FCD(Action):

    name = 'fcd'

    _questions = """
    # name of the input file
    ifile = fcd.in :: existing_file, alias=f
    """
    def run(ifile):
        theo_header.print_header('Fragment charge difference analysis')

        ifile = 'fcd.in'

        ioptions = input_options.fcd_ana_options(ifile)

        fcda = lib_diab.fcd_ana(ioptions)
        fcda.read_mos()
        fcda.read_dens()

        fcda.do_pop_ana()
        fcda.do_fcd()
