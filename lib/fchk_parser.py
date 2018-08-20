"""
Support for parsing fchk files.
A Molden file can be produced, as well, for further processing.
"""

import error_handler, file_parser, units, lib_mo, lib_struc
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

class MO_set_fchk(lib_mo.MO_set_molden):
    """
    Parse MO-related information from the fchk file.
    Not finished yet.
    """
    def read(self):
        print 'Reading MOs from fchk file %s'%self.file
        self.rfileh = open(self.file, 'r')
        while True: # loop over all lines
            try:
                line = self.rfileh.next()
            except StopIteration:
                print "Reached end of file %s"%self.rfileh.name
                break

            if 'Number of alpha electrons' in line:
                nocc = int(line.split()[-1])
            elif 'Atomic numbers' in line:
                atnos = self.fchk_list(line)
            elif 'Current cartesian coordinates' in line:
                coors = self.fchk_list(line)
                self.set_at_dicts(atnos, coors)
            elif 'Number of basis functions' in line:
                num_bas = int(line.split()[-1])

            elif 'Shell types' in line:
                shtypes = map(int, self.fchk_list(line))
            elif 'Number of primitives per shell' in line:
                noprim = map(int, self.fchk_list(line))
            elif 'Shell to atom map' in line:
                shmap = map(int, self.fchk_list(line))
            elif 'Primitive exponents' in line:
                prims = map(float, self.fchk_list(line))
            elif 'Contraction coefficients' in line:
                contr = map(float, self.fchk_list(line))
                self.set_bf_info(shtypes, noprim, shmap, prims, contr)

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
                self.occs = nocc * [2] + (self.ret_num_mo() - nocc) * [0]
            elif 'Alpha Orbital Energies' in line:
                self.ens = map(float, self.fchk_list(line))

            elif 'Beta MO coefficients' in line:
                raise error_handler.MsgError('Unrestricted calculations not supported')

    def write_molden_file(self, fname='out.mld', cfmt='% 10E', occmin=-1):
        """
        Write a file in Molden format.
        """
        self.export_AO(self.ens, self.occs, self.mo_mat.transpose(), fname, cfmt, occmin)

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

        self.header = "[Molden Format]\n"

        self.header+= "[Atoms] Angs\n"
        for iat, at in enumerate(self.at_dicts):
            atno = at['Z']
            (x, y, z) = (at['x'], at['y'], at['z'])
            self.header+= "%3s %5i %3i %12.6f %12.6f %12.6f\n"%\
                (lib_struc.Z_symbol_dict[atno], iat+1, atno, x, y, z)

    def set_bf_info(self, shtypes, noprim, shmap, prims, contr):
        num_bas={0:1, -1:3, 1:3, -2:5, 2:6, -3:7, 3:10, -4:9, 4:15}
        llab = ['S','P','D','F','G']
        self.header+= "[GTO]"

        iatl = -1
        iprim = 0
        for ish, shtype in enumerate(shtypes):
            iat = shmap[ish]
            lab = llab[abs(shtype)]
            for ibas in range(num_bas[shtype]):
                self.basis_fcts.append(lib_mo.basis_fct(iat, lab))

            if iat != iatl:
                iatl = iat
                self.header+="\n%i  0\n"%iat
            np = noprim[ish]
            self.header+= "%3s %5i 1.00\n"%(lab, np)
            for ip in range(np):
                self.header+= "%18.10E %18.10E\n"%(prims[iprim], contr[iprim])
                iprim += 1