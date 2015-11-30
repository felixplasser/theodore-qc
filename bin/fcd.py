#!/usr/bin/python
"""
Perform a fragment charge difference analysis, following
  A. A. Voityuk, N. Roesch J. Chem. Phys. 2002, 117, 5607.
"""

import theo_header, input_options, dens_ana_base, file_parser, pop_ana, error_handler
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
        rtype = self.ioptions.get('rtype')
        if rtype in ['mcscf', 'colmcscf']:
            state_list = file_parser.file_parser_col_mcscf(self.ioptions).read_fcd(self.mos)
            (self.st1, self.st2, self.trans) = state_list
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
        self.trans['pop'] = pana.ret_pop(self.trans['tden'], self.mos)
        
        if lvprt>=2:
            pop_pr = pop_ana.pop_printer(self.struc)
            pop_pr.add_pop('state 1', self.st1['pop'])
            pop_pr.add_pop('state 2', self.st2['pop'])
            pop_pr.add_pop('trans.',  self.trans['pop'])
            print pop_pr.ret_table()
            
    def do_fcd(self):
        # two entries for the two fragments
        qst1 = numpy.zeros(2)
        qst2 = numpy.zeros(2)
        qtrans = numpy.zeros(2)
        
        for atind in self.ioptions['at_lists'][0]:
            qst1[0] += self.st1['pop'][atind-1]
            qst2[0] += self.st2['pop'][atind-1]
            qtrans[0] += self.trans['pop'][atind-1]
        for atind in self.ioptions['at_lists'][1]:
            qst1[1] += self.st1['pop'][atind-1]
            qst2[1] += self.st2['pop'][atind-1]
            qtrans[1] += self.trans['pop'][atind-1]
            
        pop_pr = pop_ana.pop_printer(self.struc)
        pop_pr.add_pop('state 1', qst1)
        pop_pr.add_pop('state 2', qst2)
        pop_pr.add_pop('trans.',  qtrans)
        print pop_pr.ret_table_FCD(self.ioptions['at_lists'])


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