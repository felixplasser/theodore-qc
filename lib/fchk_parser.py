"""
Support for parsing fchk files.
Currently, this uses an auxiliary Molden file, since the Molden file has the same ordering
of the basis functions.
Later, it should create its own Molden file.
"""

import error_handler, file_parser, units, lib_mo
import numpy

class file_parser_fchk(file_parser.file_parser_base):
    """
    Read (trans)-density matrices and info from fchk file.
    The info is written with: GUI=2, state_analysis=True
    """
    def read(self, mos):
        state_list = []

        num_bas = mos.ret_num_bas()

        self.rfileh = open(self.ioptions['rfile'], 'r')
        while True: # loop over all lines
            try:
                line = self.rfileh.next()
            except StopIteration:
                print "Reached end of file %s"%self.ioptions.get('rfile')
                break

            if 'Transition DM' in line:
                print line.strip()
                state_list.append({})
                state = state_list[-1]
                words = line.split()

                state['name']    = words[0].replace('singlet','S').replace('triplet','T')
                state['exc_en']  = float(words[1]) * units.energy['eV']
                state['osc_str'] = float(words[2])
                dim = int(words[-1])

                refdim = num_bas * num_bas
                tmplist = self.fchk_list(line, refdim)
                tden_ao = numpy.reshape(map(float, tmplist), [num_bas,num_bas])
                # The tden is transformed back to the MO basis to comply with the
                #   remaining TheoDORE infrastructure
                temp = mos.CdotD(tden_ao.T, trnsp=False, inv=True)
                state['tden'] = 2**(-.5) * mos.MdotC(temp, trnsp=True, inv=True)

                # print 'DSSD   :', numpy.sum(numpy.dot(tden_ao.T,S)*numpy.dot(S,tden_ao.T))
                # print 'tden_ao:', numpy.sum(tden_ao**2)
                # print 'tden   :', numpy.sum(state['tden']**2)
                # print

            elif 'State Density' in line or 'SCF Density' in line:
                print line.strip()
                refdim = num_bas * (num_bas + 1) / 2
                tmplist = self.fchk_list(line, refdim)

                temp = numpy.zeros([num_bas, num_bas], float)
                temp[numpy.tril_indices(num_bas)] = tmplist
                sden = numpy.triu(temp.T, 1) + temp
                # print sden

                print 'DS:   ', numpy.sum(sden * mos.S)
                print

        self.rfileh.close()

        return state_list

    def fchk_list(self, line, refdim=None):
        words = line.split()
        dim = int(words[-1])
        if not refdim is None:
            assert dim == refdim, 'Inconsistent dimensions'
        tmplist = []
        while len(tmplist) < dim:
            tmplist += self.rfileh.next().split()

        return tmplist

class MO_set_fchk(lib_mo.MO_set):
    """
    Parse MO-related information from the fchk file.
    Not finished yet.
    """
    def read(self):
        self.rfileh = open(self.file, 'r')
        while True: # loop over all lines
            try:
                line = self.rfileh.next()
            except StopIteration:
                print "Reached end of file %s"%self.rfileh.name
                break

            if 'Atomic numbers' in line:
                atnos = self.fchk_list(line)
            elif 'Current cartesian coordinates' in line:
                coors = self.fchk_list(line)
                self.set_at_dicts(atnos, coors)

            elif 'Number of basis functions' in line:
                num_bas = int(line.split()[-1])

            elif 'Shell types' in line:
                shtypes = map(int, self.fchk_list(line))
            elif 'Number of primitives per shell' in line:
                noprim = self.fchk_list(line)
            elif 'Shell to atom map' in line:
                shmap = map(int, self.fchk_list(line))
                self.set_bf_info(shtypes, shmap)

            elif 'Primitive exponents' in line:
                prims = self.fchk_list(line)
            elif 'Contraction coefficients' in line:
                contr = self.fchk_list(line)
            elif 'Overlap Matrix' in line:
                refdim = num_bas * (num_bas + 1) / 2
                tmplist = self.fchk_list(line, refdim)

                temp = numpy.zeros([num_bas, num_bas], float)
                temp[numpy.tril_indices(num_bas)] = tmplist
                self.S = numpy.triu(temp.T, 1) + temp

                # print 'S'
                # print mos.S

                # temp = mos.CdotD(mos.S, trnsp=True, inv=False)
                # CSC = mos.MdotC(temp, trnsp=False, inv=False)
                # assert numpy.sum((CSC - numpy.eye(num_bas))**2) < 1.e-8, 'Inconsistent overlap matrix'

            elif 'Alpha MO coefficients' in line:
                #dim = int(line.split()[-1])
                tmplist = self.fchk_list(line)
                self.mo_mat = numpy.reshape(map(float, tmplist), [num_bas,num_bas]).T

            elif 'Beta MO coefficients' in line:
                raise error_handler.MsgError('Unrestricted calculations not supported')

    def fchk_list(self, line, refdim=None):
        words = line.split()
        dim = int(words[-1])
        if not refdim is None:
            assert dim == refdim, 'Inconsistent dimensions'
        tmplist = []
        while len(tmplist) < dim:
            tmplist += self.rfileh.next().split()

        return tmplist

    def set_at_dicts(self, atnos, coors):
        self.num_at = len(atnos)

        for iat, atno in enumerate(atnos):
            self.at_dicts.append({'Z':int(atno), 'x':float(coors[3*iat])*units.length['A'],
            'y':float(coors[3*iat+1])*units.length['A'],
            'z':float(coors[3*iat+2])*units.length['A']})

    def set_bf_info(self, shtypes, shmap):
        num_bas={0:1, -1:3, 1:3, -2:5, 2:6, -3:7, 3:10, -4:9, 4:15}
        llab = ['s','p','d','f','g']

        for ish, shtype in enumerate(shtypes):
            iat = shmap[ish]
            for ibas in range(num_bas[shtype]):
                self.basis_fcts.append(lib_mo.basis_fct(iat, llab[abs(shtype)]))