#!/usr/bin/python
"""
Driver script for transition density matrix analysis.
"""

import theo_header, lib_tden, lib_exciton, input_options, error_handler
import os, sys, time

theo_header.print_header('Transition density matrix analysis')
(tc, tt) = (time.clock(), time.time())

def ihelp():
    print " analyze_tden.py"
    print " Command line options:"
    print "  -h, -H, -help: print this help"
    print "  -ifile [dens_ana.in]: name of the input file"
    exit(0)

#--------------------------------------------------------------------------#        
# Parsing and computations
#--------------------------------------------------------------------------# 

ifile = 'dens_ana.in'

arg=sys.argv.pop(0)
while len(sys.argv)>0:
    arg = sys.argv.pop(0)
    if arg in ["-h", "-H", "-help"]:
        ihelp()
    elif arg == '-ifile':
        ifile = sys.argv.pop(0)
    else:
        raise error_handler.ElseError(arg, 'command line option')

if not os.path.exists(ifile):
    print 'Input file %s not found!'%ifile
    print 'Please create this file using theoinp or specify its location using -ifile\n'
    ihelp()
    
ioptions = input_options.tden_ana_options(ifile)

tdena = lib_tden.tden_ana(ioptions)
if 'mo_file' in ioptions: tdena.read_mos()
    
tdena.read_dens()

if'at_lists' in ioptions:
    tdena.compute_all_OmFrag()
    if ioptions['print_OmFrag']:
        tdena.fprint_OmFrag()
        
if ioptions['comp_ntos']: tdena.compute_all_NTO()

if 'RMSeh' in ioptions.get('prop_list') or 'MAeh' in ioptions.get('prop_list') or 'Eb' in ioptions.get('prop_list'):
    exca = lib_exciton.exciton_analysis()
    exca.get_distance_matrix(ioptions.get('coor_file'), ioptions.get('coor_format'))
    tdena.analyze_excitons(exca)

#--------------------------------------------------------------------------#        
# Print-out
#--------------------------------------------------------------------------# 

tdena.print_summary()

#print 'Finished at ' + time.asctime()

print "CPU time: % .1f s, wall time: %.1f s"%(time.clock() - tc, time.time() - tt)
