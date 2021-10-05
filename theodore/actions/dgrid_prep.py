#!/usr/bin/env python3

from __future__ import print_function, division
import sys

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')


class dgrid_options(input_options.write_options):
    def input(self):
        self.read_float('Mesh size (a.u.)',   'msize',   0.2)
        self.read_float('Mesh border (a.u.)', 'mborder', 4.0)
        self.read_int('Number of parallel processes', 'nproc', 1)


class DGridPrep(Action):

    name = 'dgrid_prep'

    _colt_description = 'Prepare input for DGrid'

    _user_input = """
    # Molden files
    mldfiles = :: list(existing_file)
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..input_options': 'input_options'
    })

    def run(mldfiles):

        theo_header.print_header(title=__class__._colt_description)
        print("Using the Molden files:", mldfiles)

        dopt = dgrid_options('dgrid.in')
        dopt.input()

        wfile = open('run_dgrid.bash', 'w')
        wfile.write('#!/bin/bash\n')
        wfile.write('\n###\nDGRID=dgrid\n###\n\n')

        iproc = 0
        for mldfile in mldfiles:
            # adjust the names to be compatible with dgrid
            rind   = mldfile.rfind('.')
            basen  = mldfile if rind < 0 else mldfile[:rind]
            basen2 = basen.split('/')[-1]
            bfilen = basen + '.md'
            ifilen = basen + '.inp'

            print('\nAnalyzing %s -> %s ...'%(mldfile, bfilen))

            iproc += 1
            if iproc%dopt['nproc'] == 0:
                lend = '|| exit 1'
            else:
                lend = '&'

            wfile.write('echo " *** Running %s ..."\n'%mldfile)
            wfile.write('$DGRID %s && $DGRID %s %s\n'%(mldfile, ifilen, lend))
            wfile.write('ln -s %s.md.rho_r %s.cube\n\n'%(basen, basen2))

            with open(ifilen, 'w') as ifile:
                ifile.write(f""":: dgrid_prep.py
basis={bfilen}
output=.

compute=rho
format=cube
mesh={dopt['msize']:.4f} {dopt['mborder']:.2f}\n""")

        wfile.close()
        print("File run_dgrid.bash written.\n  Run as:\n  bash run_dgrid.bash")
