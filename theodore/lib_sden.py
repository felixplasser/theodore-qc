"""
Analysis routines for state density matrices.
"""

from __future__ import print_function, division

from . import dens_ana_base, lib_mo, error_handler, pop_ana, orbkit_interface
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

        print("Number of electrons: %10.7f"%mp.sum())
        print(mp)

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
        pop_pr = pop_ana.pop_printer(self.struc)
        for dens_type in dens_types:
            pop = self.ret_general_pop(state, dens_type=dens_type)

            pop_pr.add_pop(dens_type, pop)

        print(pop_pr.ret_table())

    def print_mo_pops(self, mo_pop_type=1, lvprt=2):
        ppm = pop_ana.pop_printer_mo()
        ppm.print_mo_pops(self.mos, mo_pop_type)

    def print_all_BO(self, lvprt=2):
        """
        Print out bond order/valence information for all states.
        """
        title  = "Valence information"
        title += "\n Total valence (V_A)"
        title += "\n Free valence (F_A)"
        #title += "\n total bond order (tBO = V_A - F_A)"
        function = self.print_valence_table
        self.printer_base(title, function, lvprt)

        title  = "Bond order information"
        title += "\n <at1>-<at2> : <bond order>"
        function = self.print_BO
        self.printer_base(title, function, lvprt)

    def print_valence_table(self, state, lvprt=2, BO_data=['V_A', 'F_A']):
        pop_pr = pop_ana.pop_printer(self.struc)
        for data in BO_data:
            pop = state[data]

            pop_pr.add_pop(data, pop)

        print(pop_pr.ret_table())

    def print_BO(self, state, lvprt=2):
        for iat, BOlist in enumerate(state['BO']):
            for jat in range(iat+1, len(BOlist)):
                BOval = BOlist[jat]
                if BOval > 1.5:
                    print("{0:>3}={1:<3}:{2:>7.4f}".format(iat+1, jat+1, BOval))
                elif BOval > 0.5:
                    print("{0:>3}-{1:<3}:{2:>7.4f}".format(iat+1, jat+1, BOval))
                elif BOval > self.ioptions['min_BO']:
                    print("{0:>3}..{1:<3}:{2:>7.4f}".format(iat+1, jat+1, BOval))


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
            raise error_handler.MsgError('Population analyis type not implmented: %s'%ana_type)

        state[mp_name] = pana.ret_pop(state[dens_name], self.mos)

        return state[mp_name]

#--- Natural orbital analysis
    def compute_all_NO(self):
        """
        Analysis of natural orbitals.
        """
        if len(self.state_list) <= 1: return
        if not 'sden' in self.state_list[0]: return
        if 'nu' in self.state_list[0]: return

        jmol_orbs = self.ioptions.get('jmol_orbitals')
        if jmol_orbs:
            jmolNO = lib_mo.jmol_MOs("no")
            jmolNO.pre(ofile=self.ioptions['mo_file'])

        for state in self.state_list:
            print("NO analysis for %s"%state['name'])
            (pop, U) = self.ret_NO(state)

            if jmol_orbs:
                self.export_NOs_jmol(state, jmolNO, pop, U, minp=self.ioptions['min_occ'])

            if self.ioptions['molden_orbitals']:
                self.export_NOs_molden(state, pop, U, minp=self.ioptions['min_occ'])

        if jmol_orbs:
            jmolNO.post()

    def ret_NO(self, state):
        (pop,U) = numpy.linalg.eigh(state['sden'])

        if self.ioptions['unpaired_ana']:
            nu_v = numpy.where(pop < 1., pop, 2.-pop)
            state['nu'] = sum(nu_v)
            state['nu_den'] = numpy.dot(U,
                numpy.dot(numpy.diag(nu_v), U.T) )

            nunl_v = pop * pop * (2-pop) * (2-pop)
            state['nunl'] = sum(nunl_v)
            state['nunl_den'] = numpy.dot(U,
                numpy.dot(numpy.diag(nunl_v), U.T) )

            nel = sum(pop)
            iy0 = len(pop) - int(nel/2 + 0.5) - 1
            iy1 = len(pop) - int(nel/2 + 0.5) - 2
            state['y0'] = pop[iy0]
            state['y1'] = pop[iy1]

        return pop, U

    def export_NOs_jmol(self, state, jmolNO, pop, U, mincoeff=0.2, minp=0.01):
        Ut = U.transpose()

        jmolNO.next_set(state['name'])
        for i in range(1,len(pop)):
            pi = pop[-i]
            if pi > 2.-minp: continue
            if pi < minp: break

            jmolI = 'mo ['
            for occind in (-abs(Ut[-i])**2).argsort():
                occ = Ut[-i][occind]
                if abs(occ) < mincoeff: break
                jmolI += ' %.3f %i'%(occ,occind+1)

            jmolI += ']\n'
            jmolNO.add_mo(jmolI, "no_%s_%i"%(state['name'],i+1), pi)

    def export_NOs_molden(self, state, pop, U, minp=0.01):
        self.mos.export_MO(pop, pop, U, 'no_%s.mld'%state['name'],
                           cfmt=self.ioptions['mcfmt'], occmin=minp, alphabeta=self.ioptions['alphabeta'])

#--- Attachment / Detachment analysis

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
            print("A/D analysis for %s"%state['name'])
            (ad, W) = self.ret_NDO(state, self.state_list[0])

            if jmol_orbs:
                self.export_NDOs_jmol(state, jmolNDO, ad, W, minad=self.ioptions['min_occ'])

            if self.ioptions['molden_orbitals']:
                self.export_NDOs_molden(state, ad, W, minad=self.ioptions['min_occ'])

            if self.ioptions.get('cube_orbitals'):
                print("Calculating NDOs as cube files with orbkit.")
                oi = orbkit_interface.ok()
                oi.cube_file_creator(state, ad, W, self.mos, minlam=self.ioptions['min_occ'])

            if self.ioptions.get('pop_ana'):
                self.set_AD(state, ad, W)

        if jmol_orbs:
            jmolNDO.post()

    def compute_rho(self):
        """
        Computation of densities.
        TODO: unpaired densities.
        """
        if len(self.state_list) == 0: return
        if not 'sden' in self.state_list[0]: return

        lib_orbkit = orbkit_interface.lib_orbkit()
        cube_ids = lib_orbkit.compute_rho(self.state_list,self.mos,numproc=self.ioptions['numproc'])

    def compute_a_d_dens(self):
        # TODO: This should use NDOs
        raise error_handler.MsgError("I do not think this is correct (Felix)")
        if len(self.state_list) == 0: return
        if not 'tden' in self.state_list[0]: return
        oi = orbkit_interface.ok()
        for state in self.state_list:
            (U, lam, Vt) = self.ret_NTO(state)
            oi.compute_p_h_dens(state, U, lam, Vt, self.mos, minlam=self.ioptions['min_occ'])

    def ret_NDO(self, state, ref_state):
        """
        Compute NDOs and promotion number.
        A generalised excitation number eta is also computed here,
           cf. Barca et al. JCTC 2018, 14, 9.
        """
        dD = state['sden'] - ref_state['sden']

        (ad,W) = numpy.linalg.eigh(dD)

        state['p'] = sum(max(0., xad) for xad in ad)
        pD = sum(min(0., xad) for xad in ad)

        if abs(state['p'] + pD) > 10E-8:
            estr = 'pA + pD = %.8f != 0.'%(state['p'] + pD)
            print(' WARNING: ' + estr)

        nAA = 0.5 * numpy.sum(ref_state['sden'] * ref_state['sden'])
        nBB = 0.5 * numpy.sum(state['sden'] * state['sden'])
        nAB = 0.5 * numpy.sum(state['sden'] * ref_state['sden'])
        state["eta"] = max(nAA, nBB) - nAB

        return ad, W

    def export_NDOs_jmol(self, state, jmolNDO, ad, W, mincoeff=0.2, minad=0.05):
        Wt = W.transpose()

        jmolNDO.next_set(state['name'])
        for i, di in enumerate(ad):
            if di > -minad: break

            jmolI = 'mo color blue red\nmo ['
            for occind in (-abs(Wt[i])**2).argsort():
                occ = Wt[i][occind]
                if abs(occ) < mincoeff: break
                jmolI += ' %.3f %i'%(occ,occind+1)

            jmolI += ']\n'
            jmolNDO.add_mo(jmolI, "detach_%s_%i"%(state['name'],i+1), di)


        for i in range(1,len(ad)):
            ai = ad[-i]
            if ai < minad: break

            jmolF = 'mo color orange green\nmo ['
            for occind in (-abs(Wt[-i])**2).argsort():
                occ = Wt[-i][occind]
                if abs(occ) < mincoeff: break
                jmolF += ' %.3f %i'%(occ,occind+1)

            jmolF += ']\n'
            jmolNDO.add_mo(jmolF, "attach_%s_%i"%(state['name'],i), ai)

    def export_NDOs_molden(self, state, ad, W, minad=0.05):
        self.mos.export_MO(ad, ad, W, 'ndo_%s.mld'%state['name'],
                           cfmt=self.ioptions['mcfmt'], occmin=minad, alphabeta=self.ioptions['alphabeta'])

    def set_AD(self, state, ad, W):
        print("Computing A/D densities ...")
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

#--- Bond orders
    def compute_all_BO(self):
        """
        Compute and store the Mayer bond order matrix between the atoms.
        """
        for state in self.state_list:
            BO = self.ret_BO(state)

    def ret_BO(self, state):
        """
        Return the bond order and compute some descriptors related to the bond order.
        """
        if 'BO' in state:
            return state['BO']

        try:
            D = state['sden']
        except KeyError:
            return None

        print("Computation of the bond order matrix ...")

        temp = self.mos.CdotD(D, trnsp=False, inv=False)  # C.DAO
        DS   = self.mos.MdotC(temp, trnsp=False, inv=True) # DAO.S = C.D.C^(-1)

        # add up the contributions for the different atoms
        state['BO'] = numpy.zeros([self.mos.num_at, self.mos.num_at])

        for i in range(self.num_bas):
            iat = self.mos.basis_fcts[i].at_ind - 1
            for j in range(self.num_bas):
                jat = self.mos.basis_fcts[j].at_ind - 1
                state['BO'][iat, jat] += DS[i, j] * DS[j, i]

        QA = pop_ana.mullpop_ana().ret_pop(D, self.mos, DS)

        state['V_A'] = numpy.zeros([self.mos.num_at])
        state['F_A'] = numpy.zeros([self.mos.num_at])
        state['tBO'] = numpy.zeros([self.mos.num_at]) # direct valence
        for iat in range(self.mos.num_at):
            state['V_A'][iat] = 2 * QA[iat] - state['BO'][iat, iat]
            state['F_A'][iat] = state['V_A'][iat]
            for jat in range(self.mos.num_at):
                if jat!=iat:
                    state['F_A'][iat] -= state['BO'][iat, jat]
                    state['tBO'][iat] += state['BO'][iat, jat]

        return state['BO']
