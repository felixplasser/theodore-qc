"""
This is an interface to the external cclib Library.
http://cclib.github.io
Download and install cclib if you want to use the functions.
"""

import file_parser, error_handler, units

cclib_loaded = False
try:
    import cclib.parser
except ImportError:
    pass
else:
    cclib_loaded = True
    

class file_parser_cclib(file_parser.file_parser_base):
    def __init__(self, *args, **kwargs):
        file_parser.file_parser_base.__init__(self, *args, **kwargs)
        
    def read(self):
        state_list = []
        
        data = self.get_data(self.ioptions['rtype'], self.ioptions['rfile'])
        #attrib = data.getattributes()
        
        for ist, exc_en in enumerate(data.etenergies):
            state_list.append({})
            state = state_list[-1]
            
            state['state_ind'] = ist + 1
            state['exc_en'] = exc_en / units.energy['rcm'] * units.energy['eV']
            state['osc_str'] = data.etoscs[ist]
            state['irrep'] = data.etsyms[ist]
            
            state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
    
        return state_list
    
    def get_data(self, rtype, rfile):
        if not cclib_loaded:
            raise error_handler.MsgError('cclib not found. File type %s not supported!'%rtype)
    
        if rtype.lower() == 'gamess':
            parsec = cclib.parser.gamessparser.GAMESS
        else:
            raise error_handler.ElseError(rtype, 'cclib')
        
        return parsec(rfile).parse()
        