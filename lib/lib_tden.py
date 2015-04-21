"""
Analysis routines for transition density matrices.
"""


import dens_ana_base, Om_descriptors, lib_mo, error_handler
import numpy

numpy.set_printoptions(precision=6, suppress=True)

class tden_ana(dens_ana_base.dens_ana_base):
    """
    Analysis of transition density matrices.
    """
    # TODO: more efficient treatment for sparse matrices.
    
#--------------------------------------------------------------------------#        
# Print out
#--------------------------------------------------------------------------#     

    def print_all_tden(self, lvprt=2):
        """
        Print out all transition density matrices.
        """
        title = "Transition density matrices (MO basis)"       
        function = self.print_tden
        self.printer_base(title, function, lvprt)        
                    
    def print_tden(self, state, lvprt=2):
        tden = state['tden']
        Om = numpy.dot(tden.flatten(), tden.flatten())
    
        print "Omega = %10.7f"%Om            
        if lvprt>=2: print tden
            
#---

    def print_all_OmAt(self, lvprt=2):
        """
        Print out Omega matrices with respect to atoms.
        """
        title = "Omega matrices with respect to atoms"
        function = self.print_OmAt
        self.printer_base(title, function, lvprt)
        
    def print_OmAt(self, state, lvprt=2):
        Om, OmAt = self.ret_Om_OmAt(state)
        
        print "Omega = %10.7f"%Om            
        if lvprt>=2: print OmAt
        
#---

    def print_all_OmFrag(self, lvprt=2):
        """
        Print out Omega matrices with respect to atoms.
        """
        title = "Omega matrices with respect to fragments"
        function = self.print_OmFrag
        self.printer_base(title, function, lvprt)
        
    def print_OmFrag(self, state, lvprt=2):
        Om, OmFrag = self.ret_Om_OmFrag(state)
        
        print "Omega = %10.7f"%Om            
        if lvprt>=2: print OmFrag


    def fprint_OmFrag(self, fname="OmFrag.txt"):
        """
        Print a file containing the Omega matrix.
        """
        omf = open(fname, 'w')
        
        omf.write("%i\n"%len(self.ioptions['at_lists']))
        for state in self.state_list:
            Om, OmFrag = self.ret_Om_OmFrag(state)
            if Om == None:
                continue
            omf.write("%10s %8.5f"%(state['name'], Om))
            for el in OmFrag.flatten():
                omf.write(" %8.5f"%el)
            omf.write("\n")
            
        omf.close()
        
#---

    def print_all_Om_descriptors(self, lvprt=2, desc_list=[]):
        """
        Print Omega descriptors.
        """
        title = "Omega descriptors"
        function = self.print_Om_descriptors
        self.printer_base(title, function, lvprt, desc_list=desc_list)
        
    def print_Om_descriptors(self, state, lvprt=2, desc_list=[]):
        if not 'Om_desc' in state:
            Om, OmFrag = self.ret_Om_OmFrag(state)
            state['Om_desc'] = Om_descriptors.Om_desc_coll(Om, OmFrag)
        
        print self.ret_header_string(desc_list)
        print state['Om_desc'].ret_val_string(desc_list)
        
#---

    def print_all_exciton(self, lvprt=2):
        """
        Print out exciton analysis.
        """
        title = "Exciton analysis"
        function = self.print_exciton
        self.printer_base(title, function, lvprt)
        
    def print_exciton(self, state, lvprt=2):
        Om, OmAt = self.ret_Om_OmAt(state)
        
        print "RMS e-h sep.: %8.6f Ang"%state['RMSeh']

#--------------------------------------------------------------------------#        
# Find data
#--------------------------------------------------------------------------#     

    def ret_prop_val(self, prop, state):
        """
        Find the value of a property.
        """
        if prop in state:
            return state[prop]
        else:
            if not 'Om_desc' in state:
                Om, OmFrag = self.ret_Om_OmFrag(state)
                if Om == None: return None
                
                state['Om_desc'] = Om_descriptors.Om_desc_coll(Om, OmFrag)
            
            return state['Om_desc'].ret_desc(prop)
                
#--------------------------------------------------------------------------#        
# Operations
#--------------------------------------------------------------------------#     

    def compute_all_OmAt(self):
        """
        Computation of Omega matrices and storage in memory.
        """
        for state in self.state_list:
            Om, OmAt = self.ret_Om_OmAt(state)
            
    def ret_Om_OmAt(self, state):
        """
        Construction of the Omega matrix with respect to atoms. 
        
        formula=0: Om_mn = (DS)_mn (SD)_mn [JCTC (2012) 8, 2777]
        formula=1: Om_mn = 1/2 (DS)_mn (SD)_mn + 1/2 D_mn (SDS)_mn [JCP (2014), 141, 024106]
        formula=2: TODO - Loewdin partitioning
        """        
        if 'Om' in state and 'OmAt' in state:
            return state['Om'], state['OmAt']
        
        formula = self.ioptions.get('Om_formula')
        
        try:
            D  = state['tden']
        except KeyError:
            return None, None
        
        print "Computation of Omega matrix ..."
      # construction of intermediate matrices
        # S implicitly computed from C
        
        temp = self.mos.CdotD(D, trnsp=False, inv=False)  # C.DAO
        DS   = self.mos.MdotC(temp, trnsp=False, inv=True) # DAO.S = C.D.C^(-1)
      
        if formula==1:
            DAO = self.mos.MdotC(temp, trnsp=True, inv=False) # DAO = C.D.C^T
      
        temp = self.mos.CdotD(D, trnsp=True, inv=True)  # C^(-1,T).DAO
        SD   = self.mos.MdotC(temp, trnsp=True, inv=False)  # S.DAO = C^(-1,T).DAO.C^T
      
        if formula==1:
            # S.DAO.S = C^(-1,T).D.C^(-1)
            SDS = self.mos.MdotC(temp, trnsp=False, inv=True)

        # add up the contributions for the different atoms        
        state['Om'] = 0.
        state['OmAt'] = numpy.zeros([self.mos.num_at, self.mos.num_at])
        
        if   formula == 0:
            OmBas = DS * SD
        elif formula == 1:
            OmBas = 0.5 * (DS * SD + DAO * SDS)
        else:
            raise error_handler.MsgError("Om_formula=%i for CT numbers not implemented!"%formula)
        
        for i in xrange(self.num_bas):
            iat = self.mos.basis_fcts[i].at_ind - 1 
            for j in xrange(self.num_bas):
                jat = self.mos.basis_fcts[j].at_ind - 1
                state['Om'] += OmBas[i, j]
                state['OmAt'][iat, jat] += OmBas[i, j]
                
        return state['Om'], state['OmAt']
        
#---    

    def compute_all_OmFrag(self):
        """
        Computation of Omega matrices and storage in memory.
        """
        if not self.ioptions.has_key('at_lists'):
            print '\n WARNING: at_lists not defined - not computing CT numbers!\n'
            return
        
        for state in self.state_list:
            Om, OmFrag = self.ret_Om_OmFrag(state, self.ioptions.get('at_lists'))
            
    def ret_Om_OmFrag(self, state, at_lists=None):
        if at_lists == None:
            try:
                return state['Om'], state['OmFrag']
            except:
                return None, None
                #raise error_handler.MsgError("at_lists not specified")

        Om, OmAt = self.ret_Om_OmAt(state)
        if Om == None:
            return None, None
        
        state['OmFrag'] = numpy.zeros([len(at_lists), len(at_lists)])
        
        for A, Aatoms in enumerate(at_lists):
            for B, Batoms in enumerate(at_lists):
                for Aatom in Aatoms:
                    for Batom in Batoms:
                        state['OmFrag'][A, B] += state['OmAt'][Aatom-1, Batom-1]
                        
        return state['Om'], state['OmFrag']
#---

    def compute_all_NTO(self):
        if len(self.state_list) == 0: return
        if not 'tden' in self.state_list[0]: return
        
        jmol_orbs = self.ioptions.get('jmol_orbitals')
        if jmol_orbs:
            jmolNTO = lib_mo.jmol_MOs("nto")
            jmolNTO.pre(ofile=self.ioptions['mo_file'])
        
        for state in self.state_list:
            (U, lam, Vt) = self.ret_NTO(state)
            if jmol_orbs:
                self.export_NTOs_jmol(state, jmolNTO, U, lam, Vt)
                
            if self.ioptions['molden_orbitals']:
                self.export_NTOs_molden(state, U, lam, Vt)
            
        if jmol_orbs:
            jmolNTO.post()
            
    def ret_NTO(self, state):
        if not 'tden' in state: return None, None, None
        
        # sqrlam contains the squareroot of the singular values lambda as defined in JCP 141, 024106 (2014).
        (U, sqrlam, Vt) = numpy.linalg.svd(state['tden'])        
        lam = sqrlam * sqrlam
        
        state['PRNTO'] = lam.sum() * lam.sum() / (lam*lam).sum()
        
        return U, lam, Vt
        
    def export_NTOs_jmol(self, state, jmolNTO, U, lam, Vt, mincoeff=0.2, minlam=0.05):
        Ut = numpy.transpose(U)
        sname = state['name'].replace('(', '-').replace(')', '-')
        jmolNTO.next_set(sname)
        for i, l in enumerate(lam):
            if l < minlam: break
            
            jmolI = 'mo ['
            jmolF = 'mo ['
            
            for occind in (-abs(Ut[i])**2.).argsort():
                occ = Ut[i][occind]
                if abs(occ) < mincoeff: break
                
                jmolI += ' %.3f %i'%(occ,occind+1)
            
            for virtind in (-abs(Vt[i])**2.).argsort():
                virt = Vt[i][virtind]
                if abs(virt) < mincoeff: break
                
                jmolF += ' %.3f %i'%(virt,virtind+1)
                
            jmolI += ']\n'
            jmolNTO.add_mo(jmolI, "NTO%s_%io"%(sname,i+1), l)
            jmolF += ']\n'
            jmolNTO.add_mo(jmolF, "NTO%s_%iv"%(sname,i+1), l)
        
    def export_NTOs_molden(self, state, U, lam, Vt, mincoeff=0.2, minlam=0.01):
        """
        Export the NTOs to a molden file.
        """
        mld_name = 'nto_%s.mld'%state['name'].replace('(', '-').replace(')', '-')
        self.mos.export_NTO(lam, U, Vt, mld_name,
                           cfmt=self.ioptions['mcfmt'], occmin=minlam)
    
#---

    def analyze_excitons(self, exciton_ana):
        for state in self.state_list:
            Om, OmAt = self.ret_Om_OmAt(state)
            
            state['RMSeh'] = exciton_ana.ret_RMSeh(Om, OmAt)
            state['Eb']    = exciton_ana.ret_Eb(Om, OmAt, self.ioptions['Eb_diag'])
