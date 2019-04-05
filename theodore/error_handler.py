
from __future__ import print_function, division

class MsgError(Exception):
    def __init__(self, errmsg):
        self.errmsg = errmsg
        
    def __str__(self):
        return "\n\n  ERROR: %s"%self.errmsg

class ElseError(MsgError):
    def __init__(self, val, key):
        self.errmsg = "Option %s not implemented for '%s'!"%(val, key)
        # call e.g.
        #   raise error_handler.ElseError(ana_type, 'ana_type')

    
class NIError(MsgError):
    def __init__(self):
        self.errmsg = "Functionality no implemented yet!"
        
class PureVirtualError(MsgError):
    """
    Use this to mark pure virtual functions in an abstract base class
    These have to be redefined by the inherited class.
    -> One could use NotImplementedError instead
    """
    def __init__(self):
        self.errmsg = "Virtual function is pure!"
