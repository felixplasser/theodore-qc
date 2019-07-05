#!/usr/bin/env python3
"""
Check if a file can be read by cclib and if all the required information is available.
"""
from __future__ import print_function, division
import sys

from theodore import theo_header, cclib_interface, input_options, error_handler

theo_header.print_header('Check cclib')

print("cc_check.py <logfile> [<printlevel=1>]")
print("Check if a logfile can be parsed with cclib")

try:
    logfile = sys.argv[1]
except IndexError:
    raise error_handler.MsgError("Please enter the name of the logfile!")
if len(sys.argv) >= 3:
    lvprt = int(sys.argv[2])
else:
    lvprt = 1

ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
ioptions['rtype'] = 'cclib'
ioptions['rfile'] = logfile

ccparser = cclib_interface.file_parser_cclib(ioptions)
errcode = ccparser.check(lvprt=lvprt)

if errcode <= 1:
    print(("\n %s can be parsed by using rtype='cclib' in dens_ana.in."%logfile))
    if errcode == 0:
        print(" Conversion to Molden format also possible")
    else:
        print(" But conversion to Molden format is not possible")
else:
    print(("\n %s cannot be parsed by cclib (errcode=%i)!"%(logfile, errcode)))
