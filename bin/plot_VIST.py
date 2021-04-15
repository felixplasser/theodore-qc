#!/usr/bin/env python3
"""
Driver script for parsing NICS values and creating input for VIST plots.
"""

import sys
from theodore import theo_header, lib_NICS, error_handler, cclib_interface, input_options

def ihelp():
    print(" plot_VIST.py [options] <logfile1> <logfile2>")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -v, --vist    : VIST for only these dummy atoms, e.g. -v '0 3 5'")
    print("  -o            : Name of output file (for VMD)")
    print("  -s, --scale   : Scale factor for VIST dumb-bells")
    print("  -c, --coor    : Create coordinate files (using cclib)")
    print("  -p            : Render and plot all tensors separately")
    exit(0)

theo_header.print_header('Read NICS values and prepare VIST plot')

vlist = None
logfiles = []
do_coor = False
ofile = 'VIST.vmd'
scale = 1.
plot_all = False

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
        do_coor = True
    elif arg in ["-p"]:
        plot_all = True
    elif arg[0] == "-":
        raise error_handler.ElseError(arg, 'command line option')
    else:
        logfiles.append(arg)

if len(logfiles) == 0:
    ihelp()

open(ofile, 'w').close() # Re-initialize the file

ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
ioptions['rtype'] = 'cclib'

nv = lib_NICS.NICS_parser_g09()
for ilog, logfile in enumerate(logfiles):
    nv.read(logfile)
    nv.print_data()
    if do_coor: # create coor file to be read by VMD
        ioptions['rfile'] = logfile
        ccparser = cclib_interface.file_parser_cclib(ioptions)
        struc = cclib_interface.structure_cclib()
        struc.read_cclib(ccparser.data)
        coorf = "coor%i.xyz"%ilog
        struc.make_coord_file(file_path=coorf,file_type='Bqxyz')
        open(ofile, 'a').write("mol new %s\n"%coorf)
    nv.vmd_tensors(ofile, vlist, scale, plot_all)

# Instructions for VMD
if do_coor:
    print("""
Input file for VMD and coordinate file(s) created. Now run:
   vmd -e %s
    """%ofile)
else:
    print("""
Input file for VMD created. Now do the following:
1. Open VMD and load coordinate file
2.   File - Load Visualization State - %s
    """%ofile)