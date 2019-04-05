"""
Support for parsing fchk files.
Currently, this uses an auxiliary Molden file. Later, it should create its own Molden file.
"""

from . import error_handler, file_parser, units
import numpy

class file_parser_fchk(file_parser.file_parser_base):
    """
    Read (trans)-density matrices and info from fchk file.
    The info is written with: GUI=2, state_analysis=True
    """
    def read(self, mos):
        state_list = []
        
        rfileh = open(self.ioptions['rfile'], 'r')
        while True: # loop over all lines
            try:
                line = next(rfileh)
            except StopIteration:
                print("Reached end of file %s"%self.ioptions.get('rfile'))
                break
        
            if 'Transition density matrix' in line:
                state_list.append({})
                state = state_list[-1]
                words = line.split()
                
                state['name']    = words[0].replace('singlet','S').replace('triplet','T')
                state['exc_en']  = float(words[1]) * units.energy['eV']
                state['osc_str'] = float(words[2])
                dim = int(words[-1])
                
                num_bas = mos.ret_num_bas() # ***!!!!!!
                if dim!=num_bas*num_bas:
                    print('\n' + line)
                    print('%i != %i * %i'%(dim, num_bas, num_bas))
                    raise error_handler.MsgError('Inconsistent dimensions')
                tmplist = []
                while len(tmplist) < dim:
                    tmplist += next(rfileh).split()
                tden_ao = 2**.5 * numpy.reshape(list(map(float, tmplist)), [num_bas,num_bas])
                # The tden is transformed back to the MO basis to comply with the
                #   remaining TheoDORE infrastructure
                temp = mos.CdotD(tden_ao.T, trnsp=False, inv=True)
                state['tden'] = mos.MdotC(temp, trnsp=True, inv=True)
      
        rfileh.close()
      
        return state_list
