#!/usr/bin/env python2
"""
Use cclib and openbabel to analyze a geometry optimization.
"""

import theo_header, cclib_interface, input_options, error_handler
import openbabel
import sys

theo_header.print_header('Analyze geometry optimization')

try:
    logfile = sys.argv[1]
except IndexError:
    raise error_handler.MsgError("Please enter the name of the logfile!")

ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
ioptions['rtype'] = 'cclib'
ioptions['rfile'] = logfile

ccparser = cclib_interface.file_parser_cclib(ioptions)

obconversion = openbabel.OBConversion()
obconversion.SetOutFormat('xyz')
struc = cclib_interface.structure_cclib()
f = open('opt.xyz', 'w')

try:
    scfens = ccparser.data.scfenergies
except AttributeError:
    print ' WARNING: No SCF energies found. Quitting ...'
    sys.exit()

print 'SCF energies'
for i,scfen in enumerate(scfens):
    print '%4i: % .7f'%(i,scfen)
    
    struc.read_cclib(ccparser.data, ind=i,lvprt=0)
    lines = obconversion.WriteString(struc.mol).split('\n')
    f.write('%i\n'%ccparser.data.natom)
    f.write('Energy = % .7f\n'%scfen)
    for line in lines[2:-1]:
        f.write(line + '\n')
    
f.close()
print "Geometries written to %s"%f.name