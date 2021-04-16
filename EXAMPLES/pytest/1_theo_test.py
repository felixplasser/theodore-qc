#!/usr/bin/env python3

"""
Integration tests for TheoDORE.
Activate this by running pytest [-v/-s] in the EXAMPLES directory.
"""

# TODO: Parts of this can be moved into its own library later

import os, shutil, sys, subprocess, difflib, warnings

stddirs="pyrrole.qcadc hexatriene.colmrci fa2.ricc2 pv2p.escf pv2p.qctddft pyridine.ricc2 fa2.col fa2.rassi-libwfa fa2.terachem fa2.dftmrci fa2.cation tyrosine.ricc2-es2es biphenyl.tddftb naphth.fchk water.ricc2"
obdirs="ir_c3n3.qctddft"
cclibdirs="fa2.cclib SnH4-ecp.firefly H2S.orca"
adfdirs="fa2.adf"
faildirs="fa2.rassi"

num_run = 0
skipped = []
failed = []

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

# def test_standard():
#     run(stddirs)

def test_summary():
    """
    Print the final summary. This is printed as a warning to make it visible.
    """
    summ_str = "\n*** Integration tests finished. ***\n"
    summ_str += "Number of tests run: %i\n"%num_run
    if len(skipped) == 0:
        summ_str += "No tests skipped"
    else:
        print("Skipped tests:")
        for skip in skipped:
            summ_str += skipped + "\n"
    if len(failed) == 0:
        summ_str += "No tests failed"
    else:
        print("Failed tests:")
        for fail in failed:
            summ_str += failed + "\n"
    warnings.warn(summ_str)

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
        errname = "%s/%s"%(self.rdir, runf)
        
        ref = open(reff).readlines()
        run = open(runf).readlines()
        suffix = runf.split('.')[-1]
        if suffix == 'txt':
            for iline, line in enumerate(ref):
                if not line == run[iline]:
                    print("- ", line, end="")
                    print("+ ", run(iline), end="")
                    if not errname in failed:
                        failed.append(errname)
        else:
            diffl = list(difflib.unified_diff(ref, run, fromfile='%s (Ref.)'%rfile, tofile='%s (current)'%rfile))
            if len(diffl) > 0:
                for line in diffl:
                    print(line, end="")
                if not errname in failed:
                    failed.append(errname)
                 #assert(len(diffl) == 0)

    def run_standard(self):
        """
        Run tests in standard format.
        """
        self.prep()
        #shutil.copytree('IN_FILES', 'RUN')
        os.chdir(self.path + '/RUN')
        for ifile in sorted(os.listdir('../IN_FILES')):
            print(ifile)
            shutil.copy("../IN_FILES/"+ifile, 'dens_ana.in')
            
            finfo = ifile.split('.')
            dtype = finfo[0]
            atype = finfo[2]
            comm="analyze_%s.py"%dtype
            
            with open("analyze_%s.out"%atype, 'w') as outf:
                subprocess.run(comm, stdout=outf)
            # if dtype == 'tden':
            #     sys.argv = ['analyze_tden.py']
            #     import analyze_tden

        self.check()