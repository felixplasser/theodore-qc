#!/usr/bin/env python3

import os, subprocess
import sys
from theodore import lib_pytest
from theodore.lib_pytest import ActionFactory, mock_stdout

def test_BCyc_VIST():
    pjob = lib_pytest.pytest_job(__file__)

    with mock_stdout() as out:
        sys.argv = ['theodore', 'plot_vist', '-p', '-o', 'simple.vmd', 'neutral.log']
        ActionFactory.from_commandline()
        with open(f'run_simple.out', 'w') as fh:
            fh.write(out.read())


    with mock_stdout() as out:
        sys.argv = ['theodore', 'plot_vist', '-c', '-v', "'0 4'", 'neutral.log', 'triplet.log', '2M.log']
        ActionFactory.from_commandline()
        with open(f'run_VIST.out', 'w') as fh:
            fh.write(out.read())

    pjob.check()

print(" ".join(['theodore', 'plot_vist', '-plot_all', 'true', '-ofile', 'simple.vmd', '', 'neutral.log']))
