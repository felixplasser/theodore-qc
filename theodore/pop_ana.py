"""
Module for population analysis.
Currently only Mulliken style analysis is supported.
"""

from __future__ import print_function, division

from . import error_handler, lib_struc
import numpy

class pop_ana:
    """
    Base class for population analysis.
    """
    def ret_Deff(self, dens, mos):
        raise error_handler.PureVirtualError()

    def ret_pop(self, dens, mos, Deff=None):
        if Deff is None:
            Deff = self.ret_Deff(dens, mos)

        mp = numpy.zeros(mos.num_at)

        for ibas in range(mos.ret_num_bas()):
            iat = mos.basis_fcts[ibas].at_ind - 1
            mp[iat] += Deff[ibas, ibas]

        return mp

class mullpop_ana(pop_ana):
    """
    Mulliken population analysis.
    """
    def ret_Deff(self, dens, mos):
        """
        Compute and return the Mulliken population.
        """
        temp = mos.CdotD(dens, trnsp=False, inv=False)  # C.DAO
        DS   = mos.MdotC(temp, trnsp=False, inv=True) # DAO.S = C.D.C^(-1)

        return DS

class pop_printer:
    """
    Printer for population analysis data.
    """
    def __init__(self, struc=None):
        self.pop_types = []
        self.pops = []

        self.struc = struc

    def clear(self):
        self.pop_types = []
        self.pops = []

    ## \brief Add population data
    # \param pop_type name to be printed
    # \param pop numpy.array with data
    def add_pop(self, pop_type, pop):
        """
        Add population data to be stored in the printer class.
        """
        if pop is None: return

        self.pop_types.append(pop_type)
        self.pops.append(pop)

    def header(self, inp):
        hstr = inp
        for pop_type in self.pop_types:
            hstr += '%10s'%pop_type

        retstr  = len(hstr) * '-' + "\n"
        retstr += hstr
        retstr += "\n" + len(hstr) * '-' + "\n"

        return hstr, retstr

    def ret_table(self, labels=[]):
        """
        Return a table containing all the populations of interest.
        """
        if len(self.pop_types) == 0:
            return "  ... no population analysis data available."

        hstr, retstr = self.header('%6s'%'Atom')

        # main part
        for iat in range(len(self.pops[0])):
            if labels != []:
                retstr += '%6s'%labels[iat]
            elif self.struc is None:
                retstr += '%6i'%(iat+1)
            else:
                retstr += '%3s%3i'%(self.struc.ret_symbol(iat+1), iat+1)

            for pop in self.pops:
                retstr += '% 10.5f'%pop[iat]
            retstr += '\n'

        # sums
        retstr += len(hstr) * '-' + "\n"

        retstr += '%6s'%''
        for pop in self.pops:
            retstr += '% 10.5f'%pop.sum()

        retstr += "\n" + len(hstr) * '-' + "\n"

        return retstr

    def ret_table_Frag(self, at_lists):
        """
        Return a table over fragments.
        """
        if len(self.pop_types) == 0:
            return "  ... no population analysis data available."

        hstr, retstr = self.header('%15s'%'Fragment')

        # main part
        for i in range(len(self.pops[0])):
            if self.struc is None:
                retstr += '%15i'%(i+1)
            else:
                retstr += '%15s'%(self.struc.ret_at_list_composition(at_lists[i]))
            for pop in self.pops:
                retstr += '% 10.5f'%pop[i]
            retstr += '\n'

        # sums
        retstr += len(hstr) * '-' + "\n"

        retstr += '%15s'%''
        for pop in self.pops:
            retstr += '% 10.5f'%pop.sum()

        retstr += "\n" + len(hstr) * '-' + "\n"

        return retstr

    def ret_table_FCD(self, at_lists):
        """
        Table for FCD.
        """
        hstr, retstr = self.header('%15s'%'Fragment')

        # main part
        for i in range(len(self.pops[0])):
            if self.struc is None:
                retstr += '%15i'%(i+1)
            else:
                retstr += '%15s'%(self.struc.ret_at_list_composition(at_lists[i]))
            for pop in self.pops:
                retstr += '% 10.5f'%pop[i]
            retstr += '\n'

        # sums
        retstr += len(hstr) * '-' + "\n"

        retstr += '%15s'%'FCD'
        for pop in self.pops:
            retstr += '% 10.5f'%(pop[1]-pop[0])

        retstr += "\n" + len(hstr) * '-' + "\n"

        return retstr

class pop_printer_mo(pop_printer):
    """
    Print Mulliken populations of MOs according to atoms.
    """
    def print_mo_pops(self, mos, dosum=1, ncol=8):
        """
        Print MO populations.
        dosum: 1 - sum over atoms
               2 - sum over basis function types
        """
        if dosum==1:
            labels = []
        elif dosum==2:
            labels = mos.bf_labels
        else:
            raise error_handler.ElseError(dosum, 'dosum')

        for imo in range(mos.ret_num_mo()):
            mp = mos.ret_mo_pop(imo, dosum=dosum)
            self.add_pop('MO %i'%(imo+1), mp)

            if (imo+1) % ncol == 0:
                print(self.ret_table(labels))
                self.clear()

        if (imo+1) % ncol != 0:
            print(self.ret_table(labels))
