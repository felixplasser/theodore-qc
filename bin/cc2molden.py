#!/usr/bin/env python2
"""
Convert a log-file to Molden format with the help of cclib.
"""

import theo_header, cclib_interface, input_options, error_handler, lib_mo
import sys

theo_header.print_header('cc2molden')

print "cc2molden.py <logfile>"
print "Convert a log-file to Molden format with the help of cclib."

try:
    logfile = sys.argv[1]
except IndexError:
    raise error_handler.MsgError("Please enter the name of the logfile!")

ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
ioptions['rtype'] = 'cclib'
ioptions['rfile'] = logfile

ccparser = cclib_interface.file_parser_cclib(ioptions)

errcode = ccparser.check()
if errcode > 0:
    print " Conversion to Molden format not possible!"
    print " %s does not contain all required information"%logfile
else:
    mos = ccparser.read_mos()
    mos.write_molden_file(fname="cc.mld")

    print "\n Finished: molden format file cc.mld written."
