"""
Analysis routines for state density matrices.
"""

import dens_ana_base, lib_mo, error_handler, pop_ana
import numpy

numpy.set_printoptions(precision=6, suppress=True)

class sden_ana(dens_ana_base.dens_ana_base):
    """
    Analysis of state density matrices.
    """
    # TODO: more efficient treatment of diagonal density matrices
    # TODO: more efficient computation of Mulliken populations.

#--------------------------------------------------------------------------#        
# Print out
#--------------------------------------------------------------------------#     

    def print_all_mullpop(self, lvprt=2):
        """
        Print out all Mulliken populations.
        """
        title = "Mulliken populations"
        function = self.print_mullpop
        self.printer_base(title, function, lvprt)
        
    def print_mullpop(self, state, lvprt=2):
        mp = self.ret_general_pop(state)
        
        print "Number of electrons: %10.7f"%mp.sum()
        print mp
        
    def print_all_pop_table(self, lvprt=2):
        """
        Print out all Mulliken populations in a table.
        """
        title = "Mulliken populations"
        function = self.print_pop_table
        
        dens_types = ['state']
        if self.ioptions['unpaired_ana']: dens_types += ['nu', 'nunl']
        if self.ioptions['AD_ana']:       dens_types += ['det', 'att']
        
        self.printer_base(title, function, lvprt, dens_types=dens_types)
        
    def print_pop_table(self, state, lvprt=2, dens_types=['']):
        pop_pr = pop_ana.pop_printer()
        for dens_type in dens_types:
            pop = self.ret_general_pop(state, dens_type=dens_type)
            
            pop_pr.add_pop(dens_type, pop)
            
        print pop_pr.ret_table()
#--------------------------------------------------------------------------#        
# Operations
#--------------------------------------------------------------------------#

    def ret_general_pop(self, state, ana_type='mullpop', dens_type=''):
        """
        Return the result of a general population analysis.        
        """
        if dens_type == '' or dens_type == 'state':
            dens_name = 'sden'
            mp_name = ana_type
        else:
            dens_name = '%s_den'%dens_type
            mp_name = '%s_%s'%(ana_type, dens_type)
        
        if mp_name in state: return state[mp_name]
        if not dens_name in state: return None
        
        if ana_type == 'mullpop':
            pana = pop_ana.mullpop_ana()            
        else:
            raise error_handler.MsgError('Population analyis type no implmented: %s'%ana_type)
        
        state[mp_name] = pana.ret_pop(state[dens_name], self.mos)
            
        return state[mp_name]

#---
    
    def compute_all_AD(self):
        """
        Attachment/detachment analysis.
        """
        if len(self.state_list) <= 1: return
        if not 'sden' in self.state_list[0]: return
        
        jmol_orbs = self.ioptions.get('jmol_orbitals')
        if jmol_orbs:
            jmolNDO = lib_mo.jmol_MOs("ndo")
            jmolNDO.pre(ofile=self.ioptions['mo_file'])
        
        for state in self.state_list[1:]:
            print "A/D analysis for %s"%state['name']
            (ad, W) = self.ret_NDO(state, self.state_list[0])
            
            if jmol_orbs:
                self.export_NDOs_jmol(state, jmolNDO, ad, W)
            
            if self.ioptions['molden_orbitals']:
                self.export_NDOs_molden(state, ad, W)
            
            if self.ioptions.get('pop_ana'):    
                self.set_AD(state, ad, W)
            
        if jmol_orbs:
            jmolNDO.post()
        
    def ret_NDO(self, state, ref_state):
        dD = state['sden'] - ref_state['sden']
        
        (ad,W) = numpy.linalg.eigh(dD)
        
        state['p'] = sum(max(0., xad) for xad in ad)
        pD = sum(min(0., xad) for xad in ad)
                
        if abs(state['p'] + pD) > 10E-8:
            estr = 'pA + pD = %.8f != 0.'%(state['p'] + pD)
            print ' WARNING: ' + estr
        
        return ad, W
    
    def export_NDOs_jmol(self, state, jmolNDO, ad, W, mincoeff=0.2, minad=0.05):
        Wt = W.transpose()
        
        jmolNDO.next_set(state['name'])
        for i, di in enumerate(ad):
            if di > -minad: break
            
            jmolI = 'mo ['
            for occind in (-abs(Wt[i])**2).argsort():
                occ = Wt[i][occind]
                if abs(occ) < mincoeff: break
                jmolI += ' %.3f %i'%(occ,occind+1)
                
            jmolI += ']\n'
            jmolNDO.add_mo(jmolI, "detach_%s_%i"%(state['name'],i+1), di)
            
            
        for i in xrange(1,len(ad)):
            ai = ad[-i]
            if ai < minad: break
            
            jmolF = 'mo ['
            for occind in (-abs(Wt[-i])**2).argsort():
                occ = Wt[-i][occind]
                if abs(occ) < mincoeff: break
                jmolF += ' %.3f %i'%(occ,occind+1)
                
            jmolF += ']\n'
            jmolNDO.add_mo(jmolF, "attach_%s_%i"%(state['name'],i), ai)
            
    def export_NDOs_molden(self, state, ad, W, minad=0.05):
        self.mos.export_MO(ad, ad, W, 'ndo_%s.mld'%state['name'],
                           cfmt=self.ioptions['mcfmt'], occmin=minad)
            
    def set_AD(self, state, ad, W):
        print "Computing A/D densities ..."
        # get positive and negative indices
        #   The signs are chosen in order to make the attachment density negative and the detachment density positive
        pos = -(numpy.sign(ad)+1.)/2.
        neg =  (numpy.sign(ad)-1.)/2.
        
        state['att_den'] = numpy.dot(W,
                          numpy.dot(numpy.diag(ad*pos),
                                    numpy.transpose(W)))
        state['det_den'] = numpy.dot(W,
                          numpy.dot(numpy.diag(ad*neg),
                                    numpy.transpose(W)))
