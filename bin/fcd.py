#!/usr/bin/python
"""
Perform a fragment charge difference analysis, following
  A. A. Voityuk, N. Roesch J. Chem. Phys. 2002, 117, 5607.
"""

import theo_header, input_options, dens_ana_base, file_parser, pop_ana, error_handler, units
import sys, os, numpy

def ihelp():
    print " analyze_tden.py"
    print " Command line options:"
    print "  -h, -H, -help: print this help"
    print "  -ifile, -f [dens_ana.in]: name of the input file"
    exit(0)

class fcd_ana(dens_ana_base.dens_ana_base):
    """
    FCD analysis class.
    """
    def read_dens(self):
        self.docoup = False

        rtype = self.ioptions.get('rtype')
        if rtype in ['mcscf', 'colmcscf']:
            state_list = file_parser.file_parser_col_mcscf(self.ioptions).read_fcd(self.mos)
            (self.st1, self.st2, self.trans) = state_list
            self.docoup = True
        elif rtype in ['nos']:
            state_list = file_parser.file_parser_nos(self.ioptions).read(self.mos)
            if not len(state_list) == 2:
                raise error_handler.MsgError("Specify two states for FCD analysis")
            (self.st1, self.st2) = state_list
            self.trans = None
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
            print pop_pr.ret_table()

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
        print pop_pr.ret_table_FCD(self.ioptions['at_lists'])

        fcd1  = qst1[1] - qst1[0]
        fcd2  = qst2[1] - qst2[0]
        if self.trans: fcd12 = qtrans[1] - qtrans[0]

        Ead = (self.st2['exc_en'] - self.st1['exc_en']) * units.energy['eV']
        if self.trans:
        # full FCD formula
            etaf = 0.5 * numpy.arctan2(2*fcd12, fcd2-fcd1)
            coupf = 0.5 * Ead * numpy.sin(2*etaf)
            # explicit formula - same result
#           Ead * fcd12 / (numpy.sqrt( (fcd1 - fcd2)*(fcd1 - fcd2) + 4. * fcd12*fcd12))

        # simplified FCD formula
        fcd12_ren = (fcd2 - fcd1) / (fcd1 + fcd2)
        etas = 0.5 * numpy.arccos(fcd12_ren)
        coups = 0.5 * Ead * numpy.sin(2*etas)

        if self.trans:
            print '%15s%12s%12s'%('', 'full', 'simplified')
            print '%15s% 12.5f% 12.5f'%('Mixing angle', etaf, etas)
            if self.docoup:
                print '%15s% 12.5f% 12.5f'%('Coupling', coupf, coups)
        else:
            print '%15s%12s'%('', 'simplified')
            print '%15s% 12.5f'%('Mixing angle', etas)
            if self.docoup:
                print '%15s% 12.5f'%('Coupling', coups)


if __name__ == '__main__':
    theo_header.print_header('Fragment charge difference analysis')

    ifile = 'dens_ana.in'

    arg=sys.argv.pop(0)
    while len(sys.argv)>0:
        arg = sys.argv.pop(0)
        if arg in ["-h", "-H", "-help"]:
            ihelp()
        elif arg == '-ifile' or arg == '-f':
            ifile = sys.argv.pop(0)
        else:
            raise error_handler.ElseError(arg, 'command line option')

    if not os.path.exists(ifile):
        print 'Input file %s not found!'%ifile
        print 'Please create this file using theoinp or specify its location using -ifile\n'
        ihelp()

    ioptions = input_options.fcd_ana_options(ifile)

    fcda = fcd_ana(ioptions)
    fcda.read_mos()
    fcda.read_dens()

    fcda.do_pop_ana()
    fcda.do_fcd()