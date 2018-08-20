"""
Support for parsing fchk files.
Currently, this uses an auxiliary Molden file, since the Molden file has the same ordering
of the basis functions.
Later, it should create its own Molden file.
"""

import error_handler, file_parser, units
import numpy

class file_parser_fchk(file_parser.file_parser_base):
    """
    Read (trans)-density matrices and info from fchk file.
    The info is written with: GUI=2, state_analysis=True
    """
    def read(self, mos):
        state_list = []

        num_bas = mos.ret_num_bas()
        rfileh = open(self.ioptions['rfile'], 'r')
        while True: # loop over all lines
            try:
                line = rfileh.next()
            except StopIteration:
                print "Reached end of file %s"%self.ioptions.get('rfile')
                break

            if 'Overlap Matrix' in line:
                words = line.split()
                dim = int(words[-1])
                assert dim == num_bas * (num_bas + 1) / 2, 'Inconsistent dimensions'
                tmplist = []
                while len(tmplist) < dim:
                    tmplist += rfileh.next().split()
                temp = numpy.zeros([num_bas, num_bas], float)
                temp[numpy.tril_indices(num_bas)] = tmplist
                mos.S = numpy.triu(temp.T, 1) + temp

                # print 'S'
                # print mos.S

                temp = mos.CdotD(mos.S, trnsp=True, inv=False)
                CSC = mos.MdotC(temp, trnsp=False, inv=False)
                assert numpy.sum((CSC - numpy.eye(num_bas))**2) < 1.e-8, 'Inconsistent overlap matrix'

            if 'Transition DM' in line:
                print line.strip()
                state_list.append({})
                state = state_list[-1]
                words = line.split()

                state['name']    = words[0].replace('singlet','S').replace('triplet','T')
                state['exc_en']  = float(words[1]) * units.energy['eV']
                state['osc_str'] = float(words[2])
                dim = int(words[-1])

                if dim!=num_bas*num_bas:
                    print '\n' + line
                    print '%i != %i * %i'%(dim, num_bas, num_bas)
                    raise error_handler.MsgError('Inconsistent dimensions')
                tmplist = []
                while len(tmplist) < dim:
                    tmplist += rfileh.next().split()
                tden_ao = numpy.reshape(map(float, tmplist), [num_bas,num_bas])
                # The tden is transformed back to the MO basis to comply with the
                #   remaining TheoDORE infrastructure
                temp = mos.CdotD(tden_ao.T, trnsp=False, inv=True)
                state['tden'] = 2**(-.5) * mos.MdotC(temp, trnsp=True, inv=True)

                # print 'DSSD   :', numpy.sum(numpy.dot(tden_ao.T,S)*numpy.dot(S,tden_ao.T))
                # print 'tden_ao:', numpy.sum(tden_ao**2)
                # print 'tden   :', numpy.sum(state['tden']**2)
                # print

            if 'State Density' in line or 'SCF Density' in line:
                print line.strip()
                words = line.split()

                dim = int(words[-1])
                assert dim == num_bas * (num_bas + 1) / 2, 'Inconsistent dimensions'
                tmplist = []
                while len(tmplist) < dim:
                    tmplist += rfileh.next().split()
                temp = numpy.zeros([num_bas, num_bas], float)
                temp[numpy.tril_indices(num_bas)] = tmplist
                sden = numpy.triu(temp.T, 1) + temp
                # print sden

                print 'DS:   ', numpy.sum(sden * mos.S)
                print

        rfileh.close()

        return state_list