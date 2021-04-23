"""
General tools for using pytest.
"""

import os, shutil, subprocess, difflib
import io
import sys
from contextlib import contextmanager
from .actions import ActionFactory


class StringIO_(io.StringIO):

    def read(self):
        pos = super().tell()
        super().seek(0)
        result = super().read()
        super().seek(pos)
        return result


@contextmanager
def mock_stdout():
    old = sys.stdout
    sys.stdout = StringIO_()
    yield sys.stdout
    sys.stdout = old


class pytest_job:
    """
    Run and check job in EXAMPLES directory using pytest.
    """
    def __init__(self, cfile):
        self.wstring = ''
        self.epath = os.path.dirname(os.path.abspath(cfile))
        self.prep()

    def run_standard(self):
        """
        Run tests in standard format.
        TODO: one can include custom bash files here
        """
        os.chdir(self.epath + '/RUN')
        for ifile in sorted(os.listdir('../IN_FILES')):
            print(ifile)
            shutil.copy("../IN_FILES/"+ifile, ifile)
            col = ifile.split('.')
            dtype, atype = col[0], col[2]
            sys.argv = ['theodore', f'analyze_{dtype}',  '-ifile', ifile]
            
            with mock_stdout() as out:
                ActionFactory.from_commandline()
        assert self._check()

    def prep(self):
        """Create `RUN` dir"""
        os.chdir(self.epath)
        if os.path.exists('RUN'):
            shutil.rmtree('RUN')
        shutil.copytree('QC_FILES', 'RUN')
        os.chdir(self.epath + '/RUN')

    def _check(self):
        """ Check if there are any differences.  """
        os.chdir(self.epath + '/RUN')
        print("Checking primary output files")
        for rfile in os.listdir('../REF_FILES'):
            print("  -> " + rfile)
            self.file_diff('../REF_FILES/'+rfile, rfile)

        if len(self.wstring) > 0:
            return False
            raise pytestDiffError(self.epath + '\n\n' + self.wstring)
        return True

    def file_diff(self, reff, runf):
        """
        Check if the files are different.
         line-by-line treatment for txt files, otherwise general diff.
        """
        wstring = "\n"
        
        ref = open(reff).readlines()
        run = open(runf).readlines()
#        suffix = runf.split('.')[-1]

        # These files should match line-by-line
        if '.txt' in runf:
            for iline, line in enumerate(ref):
                if not line == run[iline]:
                    # TODO: one could add numerical thresholds here
                    self.wstring += "- " + line
                    self.wstring += "+ " + run[iline]
        # Use a general diff here
        else:
            diffl = list(difflib.unified_diff(self.diff_ignore(ref), self.diff_ignore(run), fromfile=reff, tofile=runf))
            if len(diffl) > 0:
                wstring = "\n"
                for line in diffl:
                    self.wstring += line

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


class pytestDiffError(Exception):
    def __init__(self, errmsg):
        self.errmsg = errmsg
        
    def __str__(self):
        return "\n\n *** pytest detected difference to reference results: ***\n\n %s"%self.errmsg
