import file_parser, lib_mo, error_handler, cclib_interface
import lib_mo
import numpy

class dens_ana_base:
    """
    Base class for density matrix analysis.
    """
    def __init__(self, ioptions):
        
        # state_list contains all the information about the states:
        #    quantities that are parsed from files as well as computed quantities
        self.state_list = []
        
        self.ioptions = ioptions
        
#--------------------------------------------------------------------------#        
# Input
#--------------------------------------------------------------------------#          

    def read_mos(self, lvprt=1):
        """
        Read MOs from a separate file, which is given in Molden format.
        """
        self.mos = lib_mo.MO_set_molden(file=self.ioptions.get('mo_file'))
        self.mos.read(lvprt=lvprt)
        self.read2_mos(lvprt)

    def read2_mos(self, lvprt=1):
        self.mos.compute_inverse()                    
        self.num_mo  = self.mos.ret_num_mo()
        self.num_bas = self.mos.ret_num_bas()
        
    def read_dens(self):
        """
        Read the (transition) density matrices and some supplementary information.
        """
        rtype = self.ioptions.get('rtype')
        
        if rtype=='ricc2':
            self.state_list = file_parser.file_parser_ricc2(self.ioptions).read(self.mos)
        elif rtype in ['tddft', 'escf', 'tmtddft']:
            self.state_list = file_parser.file_parser_escf(self.ioptions).read(self.mos)
        elif rtype=='libwfa':
            self.state_list = file_parser.file_parser_libwfa(self.ioptions).read()
        elif rtype=='qcadc':
            self.state_list = file_parser.file_parser_qcadc(self.ioptions).read()
        elif rtype=='qctddft':
            self.state_list = file_parser.file_parser_qctddft(self.ioptions).read(self.mos)            
        elif rtype in ['mrci', 'colmrci']:
            self.state_list = file_parser.file_parser_col_mrci(self.ioptions).read(self.mos)
        elif rtype.lower() == 'nos':
            self.state_list = file_parser.file_parser_nos(self.ioptions).read(self.mos)
        elif rtype.lower() in ['cclib', 'gamess', 'orca']:
            # these are parsed with the external cclib library
            ccli = cclib_interface.file_parser_cclib(self.ioptions)
            self.mos = ccli.read_mos()
            self.read2_mos()
            self.state_list = ccli.read(self.mos)
        else:
            raise error_handler.ElseError(rtype, 'rtype')
          
#--------------------------------------------------------------------------#          
# Output
#--------------------------------------------------------------------------#          

    def printer_base(self, title, function, lvprt=2, **kwargs):
        """
        General print-out of properties.
        """
        print
        print title
        
        for state in self.state_list:
            #print '%i%s'%(state['state_ind'],state['irrep'])            
            print state['name']
            function(state, lvprt, **kwargs)        

#---
            
    def print_summary(self):
        """
        Print a summary with information specified in prop_list.
        """
        prop_list = self.ioptions.get('prop_list')
        width = self.ioptions.get('output_prec')[0]
        ndec  = self.ioptions.get('output_prec')[1]
        oformat = '%% %i.%if'%(width, ndec)
        
        hstr  = '%-10s'%'state' + '%8s'%'dE(eV)' + '%6s'%'osc.'
        hstr +=self.ret_header_string(prop_list, width)
        
        prt_list = []
        for state in self.state_list:
            vstr  = '%-10s'%state['name'][-10:]
            vstr += ' %7.3f'%state['exc_en']
            try:
                vstr += ' %5.3f'%state['osc_str']
            except KeyError:
                vstr += '%6s'%'-'
            
            vstr += self.ret_val_string(prop_list, state, oformat)
            
            prt_list.append([state['exc_en'], vstr])
            
        prt_list.sort()
        
        ostr  = hstr + "\n"
        ostr += len(hstr) * '-' + "\n"
        for en, vstr in prt_list:
            ostr += vstr + "\n"
        
        print "\n" + ostr
        
        if 'output_file' in self.ioptions:
            ofile = self.ioptions.get('output_file')
            print "Final output copied to %s"%ofile
            open(ofile, 'w').write(ostr)
        
    def ret_header_string(self, prop_list, width=7):
        ret_str = ''
        
        for prop in prop_list:
            ret_str += '%*s'%(width, prop)
            
        return ret_str
    
    def ret_val_string(self, prop_list, state, oformat='% 7.3f'):
        ret_str = ''
        
        for prop in prop_list:
            val = self.ret_prop_val(prop, state)
                
            if val == None:
                ret_str += (len(oformat%(0.0))-1)*' ' + '-'
            else:
                try:
                    ret_str += oformat%val
                except:
                    ret_str += ' ' + val
        
        return ret_str
    
    def ret_prop_val(self, prop, state):
        """
        Find the value of a property.
        -> This is overloaded in lib_tden.py
        """
        if prop in state:
            return state[prop]
        else:
            return None

          
