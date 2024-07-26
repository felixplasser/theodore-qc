"""
Construction of Green's functions.
"""

from __future__ import print_function, division

from . import dens_ana_base, input_options, lib_mo, error_handler, orbkit_interface
import numpy

class green_options(input_options.dens_ana_options):
    """
    Input options for Green's functions.
    """
    def set_defaults(self):
        input_options.dens_ana_options.set_defaults(self)
        self['s_or_t'] = 's'
        self['minlam'] = 5.
        self['output_file']   = "green_summ.txt"

class green_ana(dens_ana_base.dens_ana_base):
    """
    Construction of Green's functions.
    """

    def compute_G(self, energies=None):
        """
        Compute Green's function.
        E - energy (default: center of HOMO-LUMO gap)
        """
        print("\n*** Computing Green's functions ... ***")

        if self.ioptions['jmol_orbitals']:
            jmolO = lib_mo.jmol_MOs("green")
            jmolO.pre(ofile=self.ioptions['mo_file'])

        for E in energies:
            # Compute the inverse Fock matrix in the Lowdin AO basis for a specific energy.
            invE = (numpy.array(self.mos.ens) - E)**(-1)
            invF = self.mos.lowdin_trans(numpy.diag(invE))
            for A, Aatoms in enumerate(self.ioptions['at_lists']):
                (U, lam, Vt) = self.ret_GNTO(invF, Aatoms)
                if self.ioptions['jmol_orbitals']:
                    self.export_NTOs_jmol({'name':'_E%.3f'%E}, jmolO, U, lam, Vt, post="_F%02i"%(A+1))

        # Post-precessing
        if self.ioptions['jmol_orbitals']:
            jmolO.post()

    # def export_jmol(self, name, jmolO, n, U, mincoeff=0.2, minn=0.05):
    #     """
    #     Export orbitals in jmol.
    #     """
    #     Ut = U.T
    #     jmolO.next_set(name)
    #     for i, ni in enumerate(n):
    #         if abs(ni) < minn:
    #             continue
    #
    #         jmolI = 'mo ['
    #         for occind in (-abs(Ut[i])**2.).argsort():
    #             occ = Ut[i][occind]
    #             if abs(occ) < mincoeff: break
    #             jmolI += ' %.3f %i'%(occ,occind+1)
    #
    #         jmolI += ']\n'
    #         jmolO.add_mo(jmolI, "G_%s_%i"%(name, i+1), ni)

    def ret_GNTO(self, invF, Aatoms):
        """
        Compute domain NTOs of the Green's function.
        """
        FA = numpy.zeros(invF.shape, float)
        for iat, ist, ien in self.mos.bf_blocks():
            if iat+1 in Aatoms:
                FA[ist:ien,:] = invF[ist:ien,:]

        (U, sqrlam, Vt) = numpy.linalg.svd(self.mos.lowdin_trans(FA, reverse=True))
        lam = sqrlam * sqrlam

        return U, lam, Vt

    def export_NTOs_jmol(self, state, jmolNTO, U, lam, Vt, mincoeff=0.2, nNTO=1, pref='G', post=''):
        Ut = numpy.transpose(U)
        sname = pref + state['name'].replace('(', '-').replace(')', '-') + post
        jmolNTO.next_set(sname)
        for i, l in enumerate(lam):
            if i >= nNTO:
                break

            jmolI = 'mo color blue red\nmo ['
            jmolF = 'mo color orange green\nmo ['

            for occind in (-abs(Ut[i])**2.).argsort():
                occ = Ut[i][occind]
                if abs(occ) < mincoeff: break

                jmolI += ' %.3f %i'%(occ,occind+1)

            for virtind in (-abs(Vt[i])**2.).argsort():
                virt = Vt[i][virtind]
                if abs(virt) < mincoeff: break

                jmolF += ' %.3f %i'%(virt,virtind+1)

            jmolI += ']\n'
            jmolNTO.add_mo(jmolI, "%s_%i_probe"%(sname,i+1), l)
            jmolF += ']\n'
            jmolNTO.add_mo(jmolF, "%s_%i_cond"%(sname,i+1), l)
