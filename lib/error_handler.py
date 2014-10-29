class MsgError(Exception):
    def __init__(self, errmsg):
        self.errmsg = errmsg
        
    def __str__(self):
        return "\n\n  ERROR: %s"%self.errmsg

class ElseError(MsgError):
    def __init__(self, option, desc):
        self.errmsg = "Option %s not implemented for %s!"%(option, desc)

    
class NIError(MsgError):
    def __init__(self):
        self.errmsg = "Functionality no implemented yet!"