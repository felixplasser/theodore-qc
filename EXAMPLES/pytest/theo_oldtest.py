#!/usr/bin/env python3

"""
Integration tests for TheoDORE.
Activate this by running pytest [-v/-s] in the EXAMPLES directory.
This is an emulation of the older script theo_test.bash
  - There are also additional testing scripts in the individual folders.
"""

# TODO: Parts of this can be moved into its own library later

from theodore import theo_header, lib_struc, lib_pytest
import os, sys, warnings

stddirs="pyrrole.qcadc hexatriene.colmrci fa2.ricc2 pv2p.escf pv2p.qctddft pyridine.ricc2 fa2.col fa2.rassi-libwfa fa2.terachem fa2.dftmrci fa2.cation tyrosine.ricc2-es2es biphenyl.tddftb naphth.fchk water.ricc2"
obdirs="ir_c3n3.qctddft"
cclibdirs="fa2.cclib SnH4-ecp.firefly H2S.orca"
adfdirs="fa2.adf"
faildirs="fa2.rassi"

tests_run = []
skipped = []

class Test2:
    def test_header(self):
        # TODO: use local path rather than THEODIR
        theo_header.print_header('TheoDORE tests')
        warnings.warn("THEODIR: " + os.environ["THEODIR"])
        sys.path.append(os.environ["THEODIR"] + "/bin")

    def run(self, rdirs):
        pjob = lib_pytest.pytest_job()
        for rdir in rdirs.split():
            tests_run.append(rdir)
            print(" *** Starting test %s ..."%rdir)
            pjob.run_standard(rdir, strict=False)
        pjob.finalise()

    def test_cclib(self):
        self.run(cclibdirs)

    def test_standard(self):
        self.run(stddirs)
    
    def test_openbabel(self):
        try:
            import openbabel
            print("obabel imported")
        except:
            warnings.warn("\n python-openbabel not found - skipping openbabel tests")
            skipped.append(obdirs)
            return
    
        self.run(obdirs)

    def test_adf(self):
        try:
            from scm.plams import KFFile
        except:
            try:
                from kf import kffile as KFFile
            except:
                warnings.warn("\n ADF not found - skipping ADF tests")
                skipped.append(adfdirs)
                return

        self.run(adfdirs)

    def test_summary(self):
        """
        Print the final summary. This is printed to pytest.out.
        """
        summ_str = "*** Integration tests finished. ***\n"
        summ_str += "Number of tests run: %i\n"%len(tests_run)
        if len(skipped) == 0:
            summ_str += "No tests skipped.\n"
        else:
            summ_str += " -> Skipped tests:\n"
            for skip in skipped:
                summ_str += skip + "\n"
        with open("%s/EXAMPLES/pytest.out"%(os.environ["THEODIR"]), 'w') as pout:
            pout.write(summ_str)