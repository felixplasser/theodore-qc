"""
This is an interface to the external cclib Library.
http://cclib.github.io
Download and install cclib if you want to use the functions.
"""

import file_parser, lib_mo, error_handler, units, lib_struc

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
        
        try:
            core_orbs = sum(self.data.coreelectrons) / 2
            if core_orbs > 0:
                print "\n WARNING: %i core orbitals detected"%core_orbs
                print " Please check the results carefully!\n"
                rect_dens = False
        except AttributeError:
            pass
            
        for ist, exc_en in enumerate(self.data.etenergies):
            state_list.append({})
            state = state_list[-1]
            
            state['state_ind'] = ist + 1
            state['exc_en'] = exc_en / units.energy['rcm'] * units.energy['eV']
            try:
                state['osc_str'] = self.data.etoscs[ist]
            except AttributeError:
                state['osc_str'] = -1.
            state['irrep'] = self.data.etsyms[ist]
            
            state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
            
            state['tden'] = self.init_den(mos, rect=rect_dens)
            #print " len:", len(state['tden']), len(state['tden'][0])
            print "\n", state['name']
            for conf in self.data.etsecs[ist]:
                [(iocc, spocc), (ivirt, spvirt), coeff] = conf
                assert(iocc < ivirt) # de-excitations cannot be handled at this point
                state['tden'][iocc, ivirt] = coeff
                
                if coeff*coeff > 0.05:
                    print " (%i->%i), coeff=% .4f"%(iocc+1, ivirt+1, coeff)
    
        return state_list
    
    def check(self, lvprt=1):
        """
        Check if the input file can be used by cclib.
        """
        errcode = 0
        
        if lvprt >= 1:
            print "\nChecking if the logfile %s can be parsed by cclib ..."%self.ioptions['rfile']
            print "\nEssential attributes:"
            
        for attr in ['mocoeffs', 'natom', 'homos', 'moenergies', 'etenergies', 'etsyms', 'etsecs']:
            chk = hasattr(self.data, attr)
            if not chk: errcode = 2
            
            if lvprt >= 1:
                print '%15s ...'%attr, chk
                
        if lvprt >= 1:
            print "\nOptional attributes:"
            
        for attr in ['etoscs', 'aooverlaps', 'mosyms']:
            chk = hasattr(self.data, attr)
            
            if lvprt >= 1:
                print '%15s ...'%attr, chk

        if lvprt >= 1:
            print "\nAttributes for creation of Molden file:"
            
        for attr in ['gbasis', 'atomcoords', 'atomnos']:
            chk = hasattr(self.data, attr)
            if not chk: errcode = max(1, errcode)
            
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
    
class MO_set_cclib(lib_mo.MO_set_molden):
    def read(self, data, lvprt=1):
        """
        Read cclib data and convert to TheoDORE format.
        """
        self.mo_mat = data.mocoeffs[0].transpose()
        if lvprt >= 1:
            print "MO-matrix read from cclib, dimension: %i x %i"%(len(self.mo_mat), len(self.mo_mat[0]))
        
        self.num_at = data.natom
            
        self.ens = data.moenergies[0]
        self.ihomo = data.homos[0]
        self.occs = (self.ihomo + 1) * [1.] + (len(self.ens) - self.ihomo - 1) * [0.]        
               
        for iat, baslist in enumerate(data.atombasis):
            for ibas in baslist:
                self.basis_fcts.append(lib_mo.basis_fct(iat, '?', '?'))
                
        try:
            self.syms = data.mosyms
        except AttributeError:
            self.syms = ['X'] * len(self.ens)
                
        try: self.S = data.aooverlaps
        except AttributeError: pass
        
        try: self.atomnos = data.atomnos
        except AttributeError: pass
        
        try: self.coords = data.atomcoords[-1]
        except AttributeError: pass
        
        try: self.gbasis = data.gbasis
        except AttributeError: pass
        
    def ret_ihomo(self):
        return self.ihomo
    
    def write_molden_file(self, fname='out.mld', cfmt='% 10E', occmin=-1):
        """
        Write a file in Molden format.
        """
        self.header = "[Molden Format]\n"

        self.header+= "[Atoms] Angs\n"
        for iat in xrange(self.num_at):
            atno = self.atomnos[iat]
            (x, y, z) = self.coords[iat]
            self.header+= "%3s %5i %3i %12.6f %12.6f %12.6f\n"%\
                (lib_struc.Z_symbol_dict[atno], iat+1, atno, x, y, z)

        self.header+= "[GTO]\n"
        for iat, at in enumerate(self.gbasis):
            self.header+= "%3i 0\n"%(iat+1)
            for bas in at:
                typ = bas[0]
                exps = bas[1]
                
                self.header+= "%3s %5i 1.00\n"%(typ, len(exps))
                
                for exp in exps:
                    self.header+= "%18.10E %18.10E\n"%exp
            self.header+= "\n"
        
        self.export_AO(self.ens, self.occs, self.mo_mat.transpose(), fname, cfmt, occmin)
