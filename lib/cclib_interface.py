"""
This is an interface to the external cclib Library.
http://cclib.github.io
Download and install cclib if you want to use the functions.
"""

import file_parser, lib_mo, error_handler, units

class file_parser_cclib(file_parser.file_parser_base):
    def __init__(self, *args, **kwargs):
        file_parser.file_parser_base.__init__(self, *args, **kwargs)
        
        self.data = self.get_data(self.ioptions['rtype'], self.ioptions['rfile'])
        
    def get_data(self, rtype, rfile):
        """
        Call cclib to read in the data.
        """
        try:
            import cclib.parser
        except ImportError:
            print "\n ERROR: Did not find the external cclib package!"
            print "Please, download and install cclib to use this functionality.\n"
            raise
               
        return cclib.parser.ccopen(rfile).parse()
    
    def read(self, mos, rect_dens=True):
        """
        Convert the information to data analyzed by TheoDORE.
        """
        state_list = []        
        
        for ist, exc_en in enumerate(self.data.etenergies):
            state_list.append({})
            state = state_list[-1]
            
            state['state_ind'] = ist + 1
            state['exc_en'] = exc_en / units.energy['rcm'] * units.energy['eV']
            state['osc_str'] = self.data.etoscs[ist]
            state['irrep'] = self.data.etsyms[ist]
            
            state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
            
            state['tden'] = self.init_den(mos, rect=rect_dens)
            for conf in self.data.etsecs[ist]:
                [(iocc, spocc), (ivirt, spvirt), coeff] = conf
                assert(iocc < ivirt) # de-excitations cannot be handled at this point
                state['tden'][iocc, ivirt] = coeff
    
        return state_list
    
    def check(self, lvprt=1):
        """
        Check if the input file can be used by cclib.
        """
        errcode = 0
        
        if lvprt >= 1:
            print "\nChecking if the logfile %s can be parsed by cclib ..."%self.ioptions['rfile']
            print "\nEssential attributes:"
            
        for attr in ['mocoeffs', 'natom', 'homos', 'moenergies', 'etenergies', 'etoscs', 'etsyms', 'etsecs']:
            chk = hasattr(self.data, attr)
            if not chk: errcode += 1
            
            if lvprt >= 1:
                print '%15s ...'%attr, chk
                
        if lvprt >= 1:
            print "\nOptional attributes:"
            
        for attr in ['aooverlaps', 'mosyms']:
            chk = hasattr(self.data, attr)
            
            if lvprt >= 1:
                print '%15s ...'%attr, chk
                
        return errcode
    
    def read_mos(self):
        """
        Return an MO_set object in TheoDORE format.
        """
        mos = MO_set_cclib(file = None)
        
        mos.read(self.data, lvprt=2)
        
        return mos
    
class MO_set_cclib(lib_mo.MO_set):
    def read(self, data, lvprt=1):
        """
        Read cclib data and convert to TheoDORE format.
        """
        self.mo_mat = data.mocoeffs[0].transpose()
        if lvprt >= 1:
            print "MO-matrix read from cclib, dimension: %i x %i"%(len(self.mo_mat), len(self.mo_mat[0]))
        
        self.num_at = data.natom
            
        self.ens = data.moenergies
        self.ihomo = data.homos[0]
        self.occs = (self.ihomo + 1) * [1.] + (len(self.ens) - self.ihomo - 1) * [0.]
                
        for iat, baslist in enumerate(data.atombasis):
            for ibas in baslist:
                self.basis_fcts.append(lib_mo.basis_fct(iat, '?', '?'))
                
        try:
            self.syms = data.mosyms
        except AttributeError:
            self.syms = ['X'] * len(self.ens)
                
        try:
            self.S = data.aooverlaps
        except AttributeError:
            pass
        
    def ret_ihomo(self):
        return self.ihomo
