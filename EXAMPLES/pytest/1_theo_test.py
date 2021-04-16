#!/usr/bin/env python3

"""
Integration tests for TheoDORE.
Activate this by running pytest [-v/-s] in the EXAMPLES directory.
"""

# TODO: Parts of this can be moved into its own library later

from theodore import theo_header, lib_struc
import os, shutil, sys, subprocess, difflib, warnings

stddirs="pyrrole.qcadc hexatriene.colmrci fa2.ricc2 pv2p.escf pv2p.qctddft pyridine.ricc2 fa2.col fa2.rassi-libwfa fa2.terachem fa2.dftmrci fa2.cation tyrosine.ricc2-es2es biphenyl.tddftb naphth.fchk water.ricc2"
obdirs="ir_c3n3.qctddft"
cclibdirs="fa2.cclib SnH4-ecp.firefly H2S.orca"
adfdirs="fa2.adf"
faildirs="fa2.rassi"

tests_run = []
skipped = []
failed = []

class Test2:
    def test_header(self):
        # TODO: use local path rather than THEODIR
        theo_header.print_header('TheoDORE tests')
        warnings.warn("THEODIR: " + os.environ["THEODIR"])
        sys.path.append(os.environ["THEODIR"] + "/bin")

    def run(self, rdirs):
        print
        for rdir in rdirs.split():
            tests_run.append(rdir)
            print(" *** Starting test %s ..."%rdir)
            tjob(rdir).run_standard()

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

        print("Avail: ", lib_struc.obabel_avail)
        assert lib_struc.obabel_avail
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
        if len(failed) == 0:
            summ_str += "No tests failed.\n"
        else:
            summ_str += " -> Failed tests:\n"
            for fail in failed:
                summ_str += fail + "\n"
        with open("%s/EXAMPLES/pytest.out"%(os.environ["THEODIR"]), 'w') as pout:
            pout.write(summ_str)

    def test_failed(self):
        """
        Raise an error messages if there were differences in any tests.
        """
        assert len(failed) == 0

###

class tjob:
    def __init__(self, rdir):
        self.rdir = rdir
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
            self.file_diff('../REF_FILES/'+rfile, rfile)

    def file_diff(self, reff, runf):
        """
        Check if the files are different.
         line-by-line treatment for txt files, otherwise general diff.
        """
        wstring = "\n"
        
        ref = open(reff).readlines()
        run = open(runf).readlines()
        suffix = runf.split('.')[-1]
        if suffix == 'txt':
            for iline, line in enumerate(ref):
                if not line == run[iline]:
                    # TODO: one could add numerical thresholds here
                    wstring += "- " + line
                    wstring += "+ " + run[iline]
        else:
            diffl = list(difflib.unified_diff(self.diff_ignore(ref), self.diff_ignore(run), fromfile=reff, tofile=runf))
            if len(diffl) > 0:
                wstring = "\n"
                for line in diffl:
                    wstring += line
                 #assert(len(diffl) == 0)

            if not wstring.strip() == '':
                failed.append("%s/%s"%(self.rdir, runf))
                warnings.warn(wstring)

    def diff_ignore(self, dlist):
        outl = []
        iglist = ["TheoDORE", "time", "rbkit", "openbabel", "capabilities"]
        for line in dlist:
            for ig in iglist:
                if ig in line:
                    break
                elif line.strip() == "":
                    break
            else:
                outl.append(line)
        return outl

    def run_standard(self):
        """
        Run tests in standard format.
        TODO: one can include custom bash files here
        """
        self.prep()
        os.chdir(self.path + '/RUN')
        for ifile in sorted(os.listdir('../IN_FILES')):
            print(ifile)
            shutil.copy("../IN_FILES/"+ifile, ifile)
            
            finfo = ifile.split('.')
            dtype = finfo[0]
            atype = finfo[2]
            comm="analyze_%s.py"%dtype
            
            with open("analyze_%s.out"%atype, 'w') as outf:
                subprocess.check_call([comm, "-f", ifile], stdout=outf)
            # if dtype == 'tden':
            #     sys.argv = ['analyze_tden.py']
            #     import analyze_tden

        self.check()
