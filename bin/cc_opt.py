#!/usr/bin/env python
"""
Use cclib and openbabel to analyze a geometry optimization.
"""

from theodore import theo_header, cclib_interface, input_options, error_handler, units
import openbabel
import sys

theo_header.print_header('Analysis of a geometry optimization')

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
    print(' WARNING: No SCF energies found. Quitting ...')
    sys.exit()

try:
    etens = ccparser.data.etenergies
    et = True
    print(' +++ Found excitation energies +++')
except AttributeError:
    etens = []
    et = False

print('\n%21s'%'SCF energies (a.u.)', end=' ')
if et: print('%15s'%'Exc. (a.u.)')
else:  print()

for i,scfen in enumerate(scfens):
    en_au = scfen/units.energy['eV']
    print('%5i:% 15.7f'%(i,en_au), end=' ')
    if et:
        try:
            en_au += etens[i]/units.energy['rcm']
            print('% 15.7f'%(en_au))
        except IndexError:
            print()
    else:
        print()

    try:
        struc.read_cclib(ccparser.data, ind=i,lvprt=0)
    except IndexError:
        pass
    else:
        lines = obconversion.WriteString(struc.mol).split('\n')
        f.write('%i\n'%ccparser.data.natom)
        f.write('Energy = % .7f\n'%(en_au))
        for line in lines[2:-1]:
            f.write(line + '\n')

f.close()
print("Geometries written to %s"%f.name)

try:
    optdone = ccparser.data.optdone
except AttributeError:
    optdone = False
if optdone:
    print("\n *** Geometry optimization converged ***")
