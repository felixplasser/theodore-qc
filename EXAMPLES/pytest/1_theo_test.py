#!/usr/bin/env python3

"""
Integration tests for TheoDORE.
Activate this by running pytest [-v/-s] in the EXAMPLES directory.
"""

# TODO: Parts of this can be moved into its own library later

import os, shutil

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

def run(rdirs):
    print
    for rdir in rdirs.split():
        print(" *** Running ", rdir)
        tjob(rdir).run_standard()

def test_cclib():
    run("fa2.cclib SnH4-ecp.firefly H2S.orca")

###

class tjob:
    def __init__(self, rdir):
        self.path = "%s/EXAMPLES/%s"%(os.environ["THEODIR"],rdir)

    def prep(self):
        os.chdir(self.path)
        if os.path.exists('RUN'):
            shutil.rmtree('RUN')
        shutil.copytree('QC_FILES', 'RUN')

    def run_standard(self):
        """
        Run tests in standard format.
        """
        self.prep()
        #shutil.copytree('IN_FILES', 'RUN')
        os.chdir(self.path + '/RUN')
        for ifile in os.listdir('../IN_FILES'):
            print(ifile)
            shutil.copy("../IN_FILES/"+ifile, '.')
            
            (dtype, dummy, atype) = ifile.split('.')
