import sys
from . import error_handler

if sys.version_info[0] != 3:
    estr  = "python3 (ideally >= v3.7) required!\n"
    estr += "  Found python version: %s"%sys.version
    raise error_handler.MsgError(estr)

from .actions import run
