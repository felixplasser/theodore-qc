#!/usr/bin/env python3

import os, subprocess
from theodore import lib_pytest

def test_BCyc_VIST():
    pjob = lib_pytest.pytest_job(__file__)

    with open("run_simple.out", 'w') as outf:
        comm = ['plot_VIST.py', '-p', '-o', 'simple.vmd', 'neutral.log']
        subprocess.check_call(comm, stdout=outf)

    with open("run_VIST.out", 'w') as outf:
        comm = ['plot_VIST.py', '-c', '-v', '0 4', 'neutral.log', 'triplet.log', '2M.log']
        subprocess.check_call(comm, stdout=outf)

    pjob.check()