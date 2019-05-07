#!/usr/bin/env python3

from __future__ import print_function, division
import os, sys
import numpy
from theodore import theo_header, lib_mo, error_handler

theo_header.print_header('Extract molden files')

print("""\
usage: Extract the hole/particle components out of a molden file
syntax: extract_molden.py <mo_file1> [<mo_file2> ...]
options: -ene          - interpret energies as occupations
         -thresh=0.001 - threshold for print-out
         -alphabeta    - use alpha/beta labels for hole/electron
""")

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

if __name__=='__main__':
    extr = extract_mld()
    mo_files = []

    args = sys.argv[1:]
    while len(args) > 0:
        arg = args.pop(0)
        if arg == '-thresh':
            extr.thresh = float(args.pop(0))
        elif arg == '-ene':
            extr.rd_ene = True
        elif arg == '-alphabeta':
            extr.decompose = False
        else:
            mo_files.append(arg)

    for mo_file in mo_files:
        extr.extract(mo_file)
