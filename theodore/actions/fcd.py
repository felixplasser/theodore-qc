#!/usr/bin/env python3
"""
Perform a fragment charge difference analysis, following
  A. A. Voityuk, N. Roesch J. Chem. Phys. 2002, 117, 5607.
"""

from __future__ import print_function, division

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    lib_diab = importer.lazy_import_as('..lib_diab', 'lib_diab')

class FCD(Action):

    name = 'fcd'

    _colt_description = 'Fragment charge difference analysis'

    _user_input = """
    # name of the input file
    ifile = fcd.in :: existing_file, alias=f
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
            '..lib_diab': 'lib_diab',
    })

    def run(ifile):
        theo_header.print_header(title=__class__._colt_description)

        ifile = 'fcd.in'

        ioptions = input_options.fcd_ana_options(ifile)

        fcda = lib_diab.fcd_ana(ioptions)
        fcda.read_mos()
        fcda.read_dens()

        fcda.do_pop_ana()
        fcda.do_fcd()
