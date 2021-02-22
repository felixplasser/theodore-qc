#!/usr/bin/env python3
"""
Driver script for parsing NICS values and creating input for VIST plots.
"""

import sys
from theodore import theo_header, lib_NICS, error_handler

def ihelp():
    print(" parse_NICS.py <logfile> [options]")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -v, --vist    : VIST for only these dummy atoms, e.g. -i '0 3 5'")
    print("  -o            : Name of output file (for VMD)")
    exit(0)

theo_header.print_header('Read NICS values and prepare VIST plot')

vlist = None
logfile = None
ofile = 'VIST.vmd'

arg = sys.argv.pop(0)
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "--help"]:
        ihelp()
    elif arg in ["-v", "--vist"]:
        arg = sys.argv.pop(0)
        vlist = [float(i) for i in arg.split()]
    elif arg in ["-o"]:
        ofile = sys.argv.pop(0)
    elif "-" in arg:
        raise error_handler.ElseError(arg, 'command line option')
    else:
        logfile = arg

if logfile is None:
    ihelp()

nv = lib_NICS.NICS_parser_g09(logfile)
nv.print_data()
nv.vmd_tensors(ofile, vlist)