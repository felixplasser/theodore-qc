#!/usr/bin/env python3

from theodore import lib_pytest

def test_standard():
    lib_pytest.pytest_job(__file__).run_standard()
