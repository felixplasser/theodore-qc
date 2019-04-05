from __future__ import print_function, division

from . import error_handler, dens_ana_base, file_parser, pop_ana, units
import numpy

"""
Routines for property based diabatization.
"""
class diab_base(dens_ana_base.dens_ana_base):
    """
    Basic routines for diabatization.
    """
    def print_summ(self, O11, O22, O12=None, Ead=None):
        """
        Print a summary.
        """
        #O_ren = (O22 - O11) / (O11 + O22)
        etas  = 0.5 * numpy.arccos((O22-O11)/2.)
        
        print("FCD analysis. Correct?")
        
        if O12==None:
            print('%15s%12s'%('', 'simplified'))
            print('%15s% 12.5f'%('Mixing angle', etas))
            if not Ead==None:
                coups = 0.5 * Ead * numpy.sin(2*etas)
                print('%15s% 12.5f'%('Coupling', coups))
        else:
            etaf = 0.5 * numpy.arctan2(2*O12, O22-O11)
            print('%15s%12s%12s'%('', 'full', 'simplified'))
            print('%15s% 12.5f% 12.5f'%('Mixing angle', etaf, etas))
            if not Ead==None:
                coups = 0.5 * Ead * numpy.sin(2*etas)
                coupf = 0.5 * Ead * numpy.sin(2*etaf)
                print('%15s% 12.5f% 12.5f'%('Coupling', coupf, coups))           

class fcd_ana(diab_base):
    """
    FCD analysis class.
    """
    def read_dens(self):
        self.Ead = None
        
        rtype = self.ioptions.get('rtype')
        if rtype in ['mcscf', 'colmcscf']:
            (self.st1, self.st2, self.trans) = file_parser.file_parser_col_mcscf(self.ioptions).read_fcd(self.mos)
            self.Ead = (self.st2['exc_en'] - self.st1['exc_en']) * units.energy['eV']
        elif rtype in ['nos']:
            state_list = file_parser.file_parser_nos(self.ioptions).read(self.mos)
            if not len(state_list) == 2:
                raise error_handler.MsgError("Specify two states for FCD analysis")
            (self.st1, self.st2) = state_list
            self.trans = None            
        elif rtype in ['mrci', 'colmrci']:
            (self.st1, self.st2, self.trans) = file_parser.file_parser_col_mrci(self.ioptions).read_fcd(self.mos)
            self.Ead = self.trans['exc_en']
        else:
            raise error_handler.ElseError(rtype, 'rtype')

        self.extra_info()

    def do_pop_ana(self, ana_type='mullpop', lvprt=2):
        if ana_type == 'mullpop':
            pana = pop_ana.mullpop_ana()
        else:
            raise error_handler.ElseError('ana_type', ana_type)

        self.st1['pop']  = pana.ret_pop(self.st1['sden'], self.mos)
        self.st2['pop']  = pana.ret_pop(self.st2['sden'], self.mos)
        if self.trans: self.trans['pop'] = pana.ret_pop(self.trans['tden'], self.mos)

        if lvprt>=2:
            pop_pr = pop_ana.pop_printer(self.struc)
            pop_pr.add_pop('state 1', self.st1['pop'])
            pop_pr.add_pop('state 2', self.st2['pop'])
            if self.trans: pop_pr.add_pop('trans.',  self.trans['pop'])
            print(pop_pr.ret_table())

    def do_fcd(self):
        # two entries for the two fragments
        qst1 = numpy.zeros(2)
        qst2 = numpy.zeros(2)
        if self.trans: qtrans = numpy.zeros(2)

        if not len(self.ioptions['at_lists']) == 2:
            raise error_handler.MsgError("Fragment definitions for two fragments required in at_lists")

        for atind in self.ioptions['at_lists'][0]:
            qst1[0] += self.st1['pop'][atind-1]
            qst2[0] += self.st2['pop'][atind-1]
            if self.trans: qtrans[0] += self.trans['pop'][atind-1]
        for atind in self.ioptions['at_lists'][1]:
            qst1[1] += self.st1['pop'][atind-1]
            qst2[1] += self.st2['pop'][atind-1]
            if self.trans: qtrans[1] += self.trans['pop'][atind-1]

        pop_pr = pop_ana.pop_printer(self.struc)
        pop_pr.add_pop('state 1', qst1)
        pop_pr.add_pop('state 2', qst2)
        if self.trans: pop_pr.add_pop('trans.',  qtrans)
        print(pop_pr.ret_table_FCD(self.ioptions['at_lists']))

        fcd1  = qst1[1] - qst1[0]
        fcd2  = qst2[1] - qst2[0]
        fcd12 = qtrans[1] - qtrans[0] if self.trans else None

        self.print_summ(fcd1, fcd2, fcd12, self.Ead)
