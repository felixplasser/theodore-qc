#!/usr/bin/env python3
"""
Use cclib and openbabel to analyze a geometry optimization.
"""

from __future__ import print_function, division
import sys

from .. import theo_header, cclib_interface, input_options, error_handler, units
import openbabel

def ihelp():
    print(" cc_opt.py <logfile>")
    print(" Command line options:")
    print("  -h, -H, --help: print this help")
    print("  -s, --scan: Analyze a relaxed scan")
    print("  -t, --thresh: Discontinuity threshold for scan")
    print("  -o, --output: Name (and path) of output file")
    exit(0)

def cc_opt():
    theo_header.print_header('Analysis of a geometry optimization or relaxed scan')
    
    
    logfile = None
    scan = False
    scan_thresh = 500
    fname = "cc_opt.xyz"
    
    arg=sys.argv.pop(0)
    while len(sys.argv)>0:
        arg = sys.argv.pop(0)
        if arg in ["-h", "-H", "--help"]:
            ihelp()
        elif arg in ["-s", "--scan"]:
            scan = True
            if fname == "cc_opt.xyz":
                fname = "cc_scan.xyz"
        elif arg in ["-t", "--thresh"]:
            scan = True
            if fname == "cc_opt.xyz":
                fname = "cc_scan.xyz"
            scan_thresh = float(sys.argv.pop(0))
        elif arg in ["-o", "--output"]:
            fname = sys.argv.pop(0)
        else:
            logfile = arg
    
    if logfile is None:
        raise error_handler.MsgError("Please enter the name of the logfile!")
    
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
    
            if abs(dE2) > scan_thresh * abs(dE1):
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
    
    f.close()
    print("Geometries written to %s"%f.name)
    
    try:
        optdone = ccparser.data.optdone
    except AttributeError:
        optdone = False
    if optdone:
        print("\n *** Geometry optimization converged ***")
