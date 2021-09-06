#!/usr/bin/env python3
"""
Use cclib and openbabel to analyze a geometry optimization.
"""

from __future__ import print_function, division
import sys

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    cclib_interface = importer.lazy_import_as('..cclib_interface', 'cclib_interface')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    units = importer.lazy_import_as('..units', 'units')




class CCOpt(Action):

    _colt_description = 'Analysis of geom. opt. or relaxed scan'

    _user_input = """
    # Logfile of quant. chemistry program
    logfile = :: existing_file
    # Analyse a relaxed scan
    scan = false :: bool, alias=s
    # Threshold for energy change (for scan)
    thresh = 500 :: float, alias=t
    # Name of output xyz file
    output  = cc_opt.xyz :: file, alias=o
    """
     
    name = 'cc_opt'

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..cclib_interface': 'cclib_interface',
            '..error_handler': 'error_handler',
            '..input_options': 'input_options',
            '..units': 'units',
    })

    def run(logfile, scan, thresh, output):
           import openbabel

           theo_header.print_header(__class__._colt_description)
           
           fname = output
           
           ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
           ioptions['rtype'] = 'cclib'
           ioptions['rfile'] = logfile
           
           ccparser = cclib_interface.file_parser_cclib(ioptions)
           
           obconversion = openbabel.OBConversion()
           obconversion.SetOutFormat('xyz')
           struc = cclib_interface.structure_cclib()
           f = open(fname, 'w')
           
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

           ngeo = 0
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
           
               # For a potential scan: find the discontinuity when one cylce is converged
               if scan:
                   try:
                       dE1 = scfen - scfens[i-1]
                   except IndexError:
                       dE1 = 1.
                   try:
                       dE2 = scfens[i+1] - scfen
                   except IndexError:
                       dE2 = 1.
           
                   if abs(dE2) > thresh * abs(dE1):
                       print('     -> Geometry written to %s (% .4f / % .4f / % .1f)'%(f.name, dE1, dE2, dE2/dE1))
                   else:
                       print('     -> Geometry skipped')
                       continue
           
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
                   ngeo += 1

           f.close()
           print("\n*** Done *** \n%i geometries written to %s"%(ngeo, f.name))
           if scan:
               print(" ... lower threshold (-t) to get more geometries.")
           
           try:
               optdone = ccparser.data.optdone
           except AttributeError:
               optdone = False
           if optdone:
               print("\n *** Geometry optimization converged ***")
