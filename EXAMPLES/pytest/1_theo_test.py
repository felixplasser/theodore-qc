#!/usr/bin/env python3

"""
Integration tests for TheoDORE.
Activate this by running pytest [-v/-s] in the EXAMPLES directory.
"""

# TODO: Parts of this can be moved into its own library later

import os, shutil, sys, subprocess, difflib

stddirs="pyrrole.qcadc hexatriene.colmrci fa2.ricc2 pv2p.escf pv2p.qctddft pyridine.ricc2 fa2.col fa2.rassi fa2.rassi-libwfa fa2.terachem fa2.dftmrci fa2.cation tyrosine.ricc2-es2es biphenyl.tddftb naphth.fchk water.ricc2"
obdirs="ir_c3n3.qctddft"
cclibdirs="fa2.cclib SnH4-ecp.firefly H2S.orca"
adfdirs="fa2.adf"

def test_header():
    from theodore import theo_header
    
    # TODO: make sure this is using the correct version here
    print
    #theo_header.print_header('TheoDORE tests')
    print("THEODIR: ", os.environ["THEODIR"])
    sys.path.append(os.environ["THEODIR"] + "/bin")

def run(rdirs):
    print
    for rdir in rdirs.split():
        print(" *** Starting test %s ..."%rdir)
        tjob(rdir).run_standard()

def test_cclib():
    run(cclibdirs)

def test_standard():
    run(stddirs)

###

class tjob:
    def __init__(self, rdir):
        self.path = "%s/EXAMPLES/%s"%(os.environ["THEODIR"],rdir)

    def prep(self):
        os.chdir(self.path)
        if os.path.exists('RUN'):
            shutil.rmtree('RUN')
        shutil.copytree('QC_FILES', 'RUN')

    def check(self):
        os.chdir(self.path + '/RUN')
        print("Checking primary output files")
        for rfile in os.listdir('../REF_FILES'):
            print("  -> " + rfile)
            ref = open('../REF_FILES/'+rfile).readlines()
            run = open(rfile).readlines()
            diffl = list(difflib.unified_diff(ref, run, fromfile='Reference', tofile='Current run'))
            for line in diffl:
                print(line, end="")
            assert(len(diffl) == 0)

    def run_standard(self):
        """
        Run tests in standard format.
        """
        self.prep()
        #shutil.copytree('IN_FILES', 'RUN')
        os.chdir(self.path + '/RUN')
        for ifile in os.listdir('../IN_FILES'):
            print(ifile)
            shutil.copy("../IN_FILES/"+ifile, 'dens_ana.in')
            
            (dtype, dummy, atype) = ifile.split('.')
            comm="analyze_%s.py"%dtype
            
            with open("analyze_%s.out"%atype, 'w') as outf:
                subprocess.run(comm, stdout=outf)
            # if dtype == 'tden':
            #     sys.argv = ['analyze_tden.py']
            #     import analyze_tden

        self.check()