#!/usr/bin/env python3

from __future__ import print_function, division
from theodore import theo_header, input_options, sys

print("""\
usage: Prepare input for DGRid
syntax: dgrid_prep.py <mo_file1> [<mo_file2> ...]
""")

class dgrid_options(input_options.write_options):
    def input(self):
        self.read_float('Mesh size (a.u.)',   'msize',   0.2)
        self.read_float('Mesh border (a.u.)', 'mborder', 4.0)
        self.read_int('Number of parallel processes', 'nproc', 1)

def run():
    mldfiles = sys.argv[1:]

    if len(mldfiles) == 0:
        print("No Molden file specified. Stopping.")
        sys.exit(0)
    else:
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
            ifile.write(""":: dgrid_prep.py
basis=%s
output=.

compute=rho
format=cube
mesh=%.4f %.2f\n"""%(bfilen, dopt['msize'], dopt['mborder']))

    wfile.close()
    print("File run_dgrid.bash written.\n  Run as:\n  bash run_dgrid.bash")

if __name__=='__main__':
    theo_header.print_header('Prepare input for DGrid')
    run()
