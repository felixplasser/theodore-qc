#!/usr/bin/env python3

from __future__ import print_function, division
import os, sys
import numpy

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_mo = importer.lazy_import_as('..lib_mo', 'lib_mo')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')


class extract_mld:
    def __init__(self, thresh=0.001, rd_ene=False, decompose=True):
        self.thresh = thresh # threshold for print-out into output file
        self.rd_ene = rd_ene # interpret energies as occupations
        self.decompose = decompose # decompose (for NDOs or NTOs)

        self.stdir = os.getcwd()

    def extract(self, mo_file):
        print("Extracting %s ..."%mo_file)
        os.chdir(self.stdir)

        mos = lib_mo.MO_set_molden(file=mo_file)
        mos.read()


        if self.decompose:
            orb_dir = "%s.dir"%mo_file
            try:
                os.chdir(orb_dir)
            except:
                os.makedirs(orb_dir)
                os.chdir(orb_dir)

            self.ex_hp(mos,  1.)
            self.ex_hp(mos, -1.)
        else:
            self.alphabeta(mos)

    def ex_hp(self, mos, sign=1.):
        """
        Extract hole or particle component
        """
        suffix = {-1.:'hole',
                   1.:'elec'}[sign]
        outfile = "%s_%s.mld"%(mos.file, suffix)

        # extract only the ones with the required sign
        mat_list = []

        ens  = []
        occs = []
        #syms = []

        if self.rd_ene:
            tmp_occs = mos.ens
        else:
            tmp_occs = mos.occs

        Ct_orig = mos.mo_mat.transpose()
        for i, occ in enumerate(tmp_occs):
            #print occ, occ*sign, self.thresh, occ*sign > self.thresh
            if occ * sign > self.thresh:
                mat_list.append(Ct_orig[i])
                ens.append(mos.ens[i])
                occs.append(occ * sign)

        Ct = numpy.array(mat_list)

        mos.export_AO(ens, occs, Ct, fname=outfile, occmin=self.thresh)
        print("  ... %s written, containing %i orbitals."%(outfile, len(ens)))

    def alphabeta(self, mos):
        """
        Use alpha/beta labels for hole/electron
        """
        outfile = "%s_ab.mld"%(mos.file)

        # extract only the ones above thresh
        mat_list = []

        ens  = []
        occs = []

        if self.rd_ene:
            tmp_occs = mos.ens
        else:
            tmp_occs = mos.occs

        Ct_orig = mos.mo_mat.transpose()
        for i, occ in enumerate(tmp_occs):
            if abs(occ) > self.thresh:
                mat_list.append(Ct_orig[i])
                ens.append(mos.ens[i])
                occs.append(occ)

        Ct = numpy.array(mat_list)

        mos.export_AO(ens, occs, Ct, fname=outfile, occmin=self.thresh, alphabeta=True)
        print("  ... %s written, containing %i orbitals."%(outfile, len(ens)))


class ExtractMolden(Action):

    _user_input = """
    # List of MO files to analyse
    mo_files = :: list(existing_file)
    # Interpret energies as occupations
    ene = False :: bool, alias=e
    # Threshold for print-out
    thresh = 0.001 :: float, alias=t
    # Use alpha/beta labels for hole/electron
    alphabeta = False :: bool, alias=ab
    """

    name = 'extract_molden'

    _colt_description = 'Extract hole/particle parts from Molden file'

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..error_handler': 'error_handler',
            '..lib_mo': 'lib_mo',
    })

    def run(mo_files, ene, thresh, alphabeta):
        theo_header.print_header(title=__class__._colt_description)

        extr = extract_mld(thresh=thresh, rd_ene=ene, decompose=not alphabeta)
        mo_files = []

        for mo_file in mo_files:
            extr.extract(mo_file)
