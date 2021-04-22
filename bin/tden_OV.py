"""
Compute the overlap between transition density matrices.
"""

from __future__ import print_function, division
from .. import theo_header, lib_tden, input_options, error_handler
import numpy
import sys, os


def ihelp():
    print(" tden_OV.py <dir1> <dir2> [<AO_OV>]")
    print(" Command line options:")
    print("  -h, -H, -help: print this help")
    print("  -ifile, -f [tden_OV.in]: name of the input file")
    print("  -ifile2, -f2 [tden_OV.in]: name of input file for the second computation")
#    print "  -m: ignore multiplicities"
    exit(0)


def tden_ov():
    ifile = 'tden_OV.in'
    ifile2 = ''
    
    largs = []
    arg=sys.argv.pop(0)
    while len(sys.argv)>0:
        arg = sys.argv.pop(0)
        if arg in ["-h", "-H", "-help"]:
            ihelp()
        elif arg == '-ifile' or arg == '-f':
            ifile = sys.argv.pop(0)
        elif arg == '-ifile2' or arg == '-f2':
            ifile2 = sys.argv.pop(0)
        elif arg[0] == '-':
            raise error_handler.ElseError(arg, 'command line option')
        else:
            largs.append(arg)
    
    if len(largs) == 2:
        dir1 = largs[0]
        dir2 = largs[1]
        AO_OV = None
    elif len(largs) == 3:
        dir1  = largs[0]
        dir2  = largs[1]
        AO_OV = largs[2]
    else:
        ihelp()
    
    if not os.path.exists(ifile):
        print('Input file %s not found!'%ifile)
        print('Please create this file using theoinp or specify its location using -ifile\n')
        ihelp()
    
    ioptions = input_options.tden_ana_options(ifile)
    theo_header.print_header('Transition density matrix overlap', ioptions=ioptions)
    
    if ifile2 == '':
        ioptions2 = ioptions
    else:
        ioptions2 = input_options.tden_ana_options(ifile2)
    
    sdir = os.getcwd()
    
    # Read info for the first job
    os.chdir(dir1)
    tdena1 = lib_tden.tden_ana(ioptions)
    if 'mo_file' in ioptions:
        tdena1.read_mos()
    tdena1.read_dens()
    os.chdir(sdir)
    
    # Read info for the second job
    os.chdir(dir2)
    tdena2 = lib_tden.tden_ana(ioptions2)
    if 'mo_file' in ioptions2:
        tdena2.read_mos()
    tdena2.read_dens()
    os.chdir(sdir)
    
    if AO_OV == None:
        print("Constructing AO-overlap matrix from MO-coefficients")
        tdena1.mos.compute_inverse()
        SMO = numpy.dot(tdena1.mos.inv_mo_mat, tdena2.mos.mo_mat)
        if ioptions['lvprt'] >= 2:
            SAO = numpy.dot(tdena1.mos.inv_mo_mat.T, tdena1.mos.inv_mo_mat)
    else:
        SAO  = numpy.array([[float(s) for s in line.split()] for line in open(AO_OV, 'r').readlines()[1:]])
    
        CS  = numpy.dot(tdena1.mos.mo_mat.T, SAO)
        SMO = numpy.dot(CS, tdena2.mos.mo_mat)
    
    if ioptions['lvprt'] >= 2:
        print("AO-overlap matrix:", SAO.shape)
        print(SAO)
        print()
    
        print("MO-overlap matrix:", SMO.shape)
        print(SMO)
        print()
    
    print("        ", end=' ')
    for state2 in tdena2.state_list:
        print("   |%7s>"%state2['name'], end=' ')
    print()
    
    for state1 in tdena1.state_list:
        print("<%7s|"%state1['name'], end=' ')
        tden1 = state1['tden']
        DS1 = numpy.dot(tden1, SMO)
        #print DS1.shape
        mdim = tden1.shape[0]
        SDS1 = numpy.dot(SMO[:mdim,:mdim], DS1)
        #print SDS1.shape
        for state2 in tdena2.state_list:
            #if (not 'mult' in state1 or not 'mult' in state2) or (state1['mult'] == state2['mult']):
                tden2 = state2['tden']
                #print tden2.shape
                OV = sum((SDS1 * tden2).flatten())
                print(" % .8f"%OV, end=' ')
            #else:
                #print " % .1f       "%0.,
        print()
