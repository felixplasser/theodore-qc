"""
Compute the overlap between transition density matrices.
"""

from __future__ import print_function, division
from .actions import Action
import numpy
import os
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_tden = importer.lazy_import_as('..lib_tden', 'lib_tden')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')


class TDenOv(Action):

    name = 'tden_ov'

    _colt_description = 'Transition density matrix overlap'

    _user_input = """
    dir1 = :: existing_folder
    dir2 = :: existing_folder
    ao_ov = :: existing_file, optional
    # name of the input file
    ifile = tden_OV.in :: existing_file, alias=f
    # name of input file for the second computation
    ifile2 = :: existing_file, optional, alias=f2
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_tden': 'lib_tden',
            '..error_handler': 'error_handler',
            '..input_options': 'input_options',
    })

    def run(dir1, dir2, ao_ov, ifile, ifile2):
        ifile = 'tden_OV.in'
        ifile2 = ''
        
        ioptions = input_options.tden_ana_options(ifile)
        theo_header.print_header(title=__class__._colt_description, ioptions=ioptions)
        
        if ifile2 is None:
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
