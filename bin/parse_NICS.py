#!/usr/bin/env python3
"""
Driver script for parsing NICS values and creating input for VIST plots.
"""

import sys
from theodore import theo_header, lib_NICS, error_handler

def ihelp():
    print(" parse_NICS.py [options] <logfile1> <logfile2>")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -v, --vist    : VIST for only these dummy atoms, e.g. -i '0 3 5'")
    print("  -o            : Name of output file (for VMD)")
    print("  -s, --scale   : Scale factor VIST dumb-bells")
    print("  -c, --coor    : Coordinate file associated to logfile (should be used several times for several logfiles)")
    exit(0)

theo_header.print_header('Read NICS values and prepare VIST plot')

vlist = None
logfiles = []
coorfiles = []
ofile = 'VIST.vmd'
scale = 1.

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
    elif arg in ["-s", "--scale"]:
        scale = float(sys.argv.pop(0))
    elif arg in ["-c", "--coor"]:
        coorfiles.append(sys.argv.pop(0))
    elif arg[0] == "-":
        raise error_handler.ElseError(arg, 'command line option')
    else:
        logfiles.append(arg)

if len(logfiles) == 0:
    ihelp()

open(ofile, 'w').close() # Re-initialize the file

nv = lib_NICS.NICS_parser_g09()
for ilog, logfile in enumerate(logfiles):
    nv.read(logfile)
    nv.print_data()
    try:
        open(ofile, 'a').write("mol new %s\n"%coorfiles[ilog])
    except IndexError:
        pass
    nv.vmd_tensors(ofile, vlist, scale)
