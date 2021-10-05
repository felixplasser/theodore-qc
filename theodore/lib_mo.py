"""
Handling and manipulation of MO-coefficients.
"""

from __future__ import print_function, division

from . import error_handler, lib_file, units
import numpy

class MO_set:
    """
    Main class that contains orbital information.
    """
    def __init__(self, file, read=False):
        self.file = file

        self.header = '' # header info for print-out
        self.num_at = 0
        self.basis_fcts = [] # info about basis functions
        self.bf_labels = [] # list of basis function labels
        self.at_dicts = [] # info about atoms: {'Z':, 'x':, 'y':, 'z':}

        self.S = None
        self.mo_mat = None
        self.inv_mo_mat = None
        self.lowdin_mat = None
        self.Sinv2 = None # S^(-1/2)

        if read:
            self.read()

    def read(self, *args, **kwargs):
        """
        Read MOs from external file.
        """
        raise error_handler.PureVirtualError()

    def comp_inv_lowdin(self, Om_formula, lvprt=1):
        """
        Compute either the inverse or Lowdin matrix depending on Om_formula.
        """
        if Om_formula <= 1:
            self.compute_inverse(lvprt)
        elif Om_formula == 2:
            self.compute_lowdin_mat(lvprt)

    def compute_inverse(self, lvprt=1):
        """
        Compute the inverse of the MO matrix.
        If the matrix is not square, the pseudoinverse is computed.
        """

        # Preferably, the overlap matrix should be used to avoid explicit inversion
        if not self.S is None:
            if lvprt >= 1:
                print(" ... inverse computed as: C^T.S")
            self.inv_mo_mat = numpy.dot(self.mo_mat.T, self.S)
        elif len(self.mo_mat) == len(self.mo_mat[0]):
            if lvprt >= 1:
                print(" ... inverting C")
            try:
                self.inv_mo_mat = numpy.linalg.inv(self.mo_mat)
            except:
                if lvprt >= 1:
                    print(" WARNING: inversion failed.")
                    print('  Using the Moore-Penrose pseudo inverse instead.')
                self.inv_mo_mat = numpy.linalg.pinv(self.mo_mat)
        else:
            if lvprt >= 1:
                print('MO-matrix not square: %i x %i'%(len(self.mo_mat),len(self.mo_mat[0])))
                print('  Using the Moore-Penrose pseudo inverse.')
            self.inv_mo_mat = numpy.linalg.pinv(self.mo_mat)

    def compute_lowdin_mat(self, lvprt=1):
        """
        Compute the transformation matrix for Lowdin orthogonalization.
        """
        print("Performing Lowdin orthogonalization")

        (U, sqrlam, Vt) = numpy.linalg.svd(self.mo_mat)

        if Vt.shape[0] == U.shape[1]:
            self.lowdin_mat = numpy.dot(U, Vt)
            # The matrix S^(-1/2) for backtransformation
            self.Sinv2 = numpy.dot(U*sqrlam, U.T)
        elif Vt.shape[0] < U.shape[1]:
            Vts = Vt.shape[0]
            print('  MO-matrix not square: %i x %i'%(len(self.mo_mat),len(self.mo_mat[0])))
            self.lowdin_mat = numpy.dot(U[:,:Vts], Vt)
            self.Sinv2 = numpy.dot(U[:,:Vts]*sqrlam, U.T[:Vts,:])
        else:
            raise error_handler.ElseError('>', 'Lowdin ortho')

    def ret_mo_mat(self, trnsp=False, inv=False):
        """
        Return the MO matrix, possibly transposed and/or inverted.
        """
        if inv and self.inv_mo_mat is None:
            self.compute_inverse()

        if not trnsp and not inv:
            return self.mo_mat
        elif trnsp and not inv:
            return self.mo_mat.transpose()
        elif not trnsp and inv:
            return self.inv_mo_mat
        elif trnsp and inv:
            return self.inv_mo_mat.transpose()

    def ret_mo_pop(self, imo, dosum=0):
        """
        Return a vector for Mulliken population analysis of the MO.
        dosum: 0 - return the vector directly
               1 - sum over atoms
               2 - sum over basis function types
        """
        MOi    = self.mo_mat[:,imo]
        invMOi = self.inv_mo_mat[imo,:]

        CCinv = MOi * invMOi

        if dosum==0:
            return CCinv

        elif dosum==1:
            mp = numpy.zeros(self.num_at)
            for ibas in range(self.ret_num_bas()):
                iat = self.basis_fcts[ibas].at_ind - 1
                mp[iat] += CCinv[ibas]
            return mp

        elif dosum==2:
            mp = numpy.zeros(len(self.bf_labels))
            for ibas in range(self.ret_num_bas()):
                ilab = self.bf_labels.index(self.basis_fcts[ibas].label())
                mp[ilab] += CCinv[ibas]
            return mp

    def ret_ihomo(self):
        """
        Return the HOMO index (starting with 0).
        This only works if no fractional occupation numbers are present and only
           for spin-restricted orbitals.
        """
        return self.occs.index(0.) - 1

    def ret_FMO_energy(self, occ=True, offset=0):
        """
        Return the energy of frontier MOs (no assumption about ordering).
        To get the HOMO-2: set occ=True, offset=2
        """
        self.sort_ens = []
        if occ:
            for ien, en in enumerate(self.ens):
                if self.occs[ien] >= 1:
                    self.sort_ens.append(en)
            self.sort_ens.sort()
            return self.sort_ens[-offset-1]
        else:
            for ien, en in enumerate(self.ens):
                if self.occs[ien] < 1.:
                    self.sort_ens.append(en)
            self.sort_ens.sort()
            return self.sort_ens[offset]

    def ret_num_mo(self):
        return len(self.mo_mat[0])

    def ret_num_bas(self):
        return len(self.mo_mat)

    def ret_eo(self, imo):
        return self.ens[imo], self.occs[imo]

    def ret_sym(self, imo):
        try:
            return self.syms[imo]
        except:
            print("\nNo entry for imo=%i"%imo)
            raise

    def set_ens_occs(self):
        """
        In the case of Q-Chem the occupations are actually written into the energy field.
        This subroutine can be used to switch this.
        """
        self.occs = self.ens

    def MdotC(self, M, trnsp=True, inv=False):
        """
        Right-multiplication of matrix M with the MO-coefficients.
        """
        try:
            return numpy.dot(M, self.ret_mo_mat(trnsp, inv))
        except:
            print("M: %i x %i"%(len(M), len(M[0])))
            print("C: %i x %i"%(len(self.ret_mo_mat(trnsp, inv)), len(self.ret_mo_mat(trnsp, inv)[0])))
            raise

    def CdotD(self, D, trnsp=False, inv=False):
        """
        Left-multiplication of matrix D with the MO-coefficients.
        Optionally, D can be a rectangular matrix of dimension occ x (occ + virt).
        """
        if self.ret_num_mo() == len(D):
            return numpy.dot(self.ret_mo_mat(trnsp, inv), D)

        # Handling of special cases
        elif self.ret_num_mo() > len(D):
            if not trnsp and not inv:
                # take only the occ. subblock
                Csub = self.mo_mat[:,:len(D)]
                return numpy.dot(Csub, D)
            elif trnsp and inv:
                # take only the occ. subblock
                Csub = self.inv_mo_mat[:len(D)]
                return numpy.dot(Csub.transpose(), D)
            else:
                raise error_handler.ElseError('"transpose xor inverse"', 'CdotD')

        elif self.ret_num_mo() < len(D):
                print("\n WARNING: C/D mismatch")
                print(" C: %i x %i"%(self.ret_num_mo(), self.ret_num_bas()))
                print(" D: %i x %i"%(len(D), len(D[0])))
#                raise error_handler.ElseError('C < D', 'MO matrix')

                Dsub = D[:self.ret_num_mo()]
                return numpy.dot(self.ret_mo_mat(trnsp, inv), Dsub)

    def OmBas_Mulliken(self, D, formula):
        """
        Compute OmBas in a Mulliken-style analysis.
        """
        # construction of intermediate matrices
           # S implicitly computed from C
        tempA = self.CdotD(D, trnsp=False, inv=False)  # C.D
        DS   = self.MdotC(tempA, trnsp=False, inv=True) # DAO.S = C.D.C^(-1)
        tempB = self.CdotD(D, trnsp=True, inv=True)  # C^(-1,T).D
        SD   = self.MdotC(tempB, trnsp=True, inv=False)  # S.DAO = C^(-1,T).DAO.C^T

        if   formula == 0:
            return DS * SD
        elif formula == 1:
            DAO = self.MdotC(tempA, trnsp=True, inv=False) # DAO = C.D.C^T
            # S.DAO.S = C^(-1,T).D.C^(-1)
            SDS = self.MdotC(tempB, trnsp=False, inv=True)
            return 0.5 * (DS * SD + DAO * SDS)
        else:
            raise error_handler.ElseError(formula, 'Om formula')

    def comp_OmAt(self, OmBas):
        """
        Compute the Omega matrix wrt atoms.
        Use an optimised algorithm via bf_blocks.
        """
        OmAt = numpy.zeros([self.num_at, self.num_at])
        bf_blocks = self.bf_blocks()
        for iat, ist, ien in bf_blocks:
            for jat, jst, jen in bf_blocks:
                OmAt[iat, jat] = numpy.sum(OmBas[ist:ien, jst:jen])

        return OmAt

    def lowdin_trans(self, D, reverse=False):
        """
        MO-AO transformation and Lowdin orthogonalization by using
           S^0.5 C = U V^T
        """
        if not reverse:
            Lmat = self.lowdin_mat
            Rmat = self.lowdin_mat.T
        else:
            Lmat = self.lowdin_mat.T
            Rmat = self.lowdin_mat

        DUTT = numpy.dot(D, Rmat)
        if Lmat.shape[1] == DUTT.shape[0]:
            return numpy.dot(Lmat, DUTT)
        elif Lmat.shape[1] > DUTT.shape[0]:
            return numpy.dot(Lmat[:,:DUTT.shape[0]], DUTT)
        else:
            raise error_handler.ElseError('<', 'Lowdin trans')

    def lowdin_AO_trans(self, D):
        """
        Transformation from Lowdin basis to normal AO basis using
          DAO = S^(-0.5) DLow S^(-0.5) and
          S^(-0.5) = C (U V^T)^T
        """
        #if self.Sinv2 is None:
        #    self.Sinv2 = self.CdotD(self.lowdin_mat.T)

        return numpy.dot(self.Sinv2, numpy.dot(D, self.Sinv2))

    def export_MO(self, ens, occs, U, *args, **kwargs):
        """
        Exports NO, NDO etc. coefficients given in the MO basis.
        """
        mo_mat = self.CdotD(U, trnsp=False, inv=False)

        self.export_AO(ens, occs, mo_mat.transpose(), *args, **kwargs)

    def export_NTO(self, lam, U, Vt, *args, **kwargs):
        """
        Exports NTO coefficients given in the MO basis.
        """
        U_mat_t = self.CdotD(U,  trnsp=False, inv=False).transpose()
        V_mat_t = self.MdotC(Vt, trnsp=True,  inv=False)

        UV_t = [iU for iU in reversed(U_mat_t)] + [iV for iV in V_mat_t]

        lam2 = [-vlam for vlam in reversed(lam)] + [vlam for vlam in lam] + [0.] * (len(Vt) - len(lam))

        self.export_AO(lam2, lam2, UV_t, *args, **kwargs)

    def export_AO(self, *args, **kwargs):
        raise error_handler.PureVirtualError()

    def bf_blocks(self, bf_list = None, num_bas = None):
        """
        Return a list with the start and end indices for basis functions on the different atoms.
        [(iat, ist, ien), ...]
        """
        if bf_list is None:
            bf_list = self.basis_fcts
        if num_bas is None:
            num_bas = self.ret_num_bas()

        bf_blocks = []
        iat_old = 0
        ist = 0
        for ibas in range(num_bas):
            iat = bf_list[ibas].at_ind - 1
            if iat != iat_old:
                bf_blocks.append((iat_old, ist, ibas))
                iat_old = iat
                ist = ibas
        bf_blocks.append((iat, ist, num_bas))

        return bf_blocks

    def symsort(self, irrep_labels, Om_formula, sepov=True):
        """
        Sort MOs by symmetry (in case they are sorted by energy).
        This is more of a hack than a clean and stable routine...
        """
        print("Sorting MOs by symmetry")

        occorbs = {}
        virtorbs = {}
        for il in irrep_labels:
            occorbs[il]  = []
            virtorbs[il] = []

        for imo, sym in enumerate(self.syms):
            for il in irrep_labels:
                if il in sym:
                    if self.occs[imo] == 0.:
                        virtorbs[il].append((sym, imo))
                    else:
                        occorbs[il].append((sym, imo))

        T = numpy.zeros([self.ret_num_mo(), self.ret_num_mo()], int)

        orblist = []
        if sepov:
            for il in irrep_labels:
                orblist += occorbs[il]
            for il in irrep_labels:
                orblist += virtorbs[il]
        else:
            raise error_handler.ElseError('False', 'sepov')

        self.syms = []
        for jmo, orb in enumerate(orblist):
            T[jmo, orb[1]] = 1
            self.syms.append(orb[0])

        assert(jmo==self.ret_num_mo()-1)

        self.mo_mat = numpy.dot(self.mo_mat, T.transpose())
        self.comp_inv_lowdin(Om_formula)

class MO_set_molden(MO_set):
    def export_AO(self, ens, occs, Ct, fname='out.mld', cfmt='% 10E', occmin=-1, alphabeta=False):
        """
        Export coefficients given already in the AO basis to molden file.

        Ct can either be a list or numpy array with the coefficients.
        """
        mld = open(fname, 'w')
        mld.write(self.header)
        mld.write('[MO]\n')

        for imo in range(len(Ct)):
            if abs(occs[imo]) < occmin: continue

            mld.write(' Sym= X\n')
            mld.write(' Ene= %f\n'%ens[imo])
            if alphabeta:
                if occs[imo] < 0:
                    mld.write(' Spin= Alpha\n')
                    mld.write(' Occup= %f\n'%-occs[imo])
                else:
                    mld.write(' Spin= Beta\n')
                    mld.write(' Occup= %f\n'%occs[imo])
            else:
                mld.write(' Spin= Alpha\n')
                mld.write(' Occup= %f\n'%occs[imo])
            for ibf, coeff in enumerate(Ct[imo]):
                fmtstr = '%10i   '+cfmt+'\n'
                mld.write(fmtstr%(ibf+1, coeff))

        mld.close()

    def ret_coeffs(self, occmin=-1., occmax=100., eneocc=False, sym='X'):
        outstr = ''

        for imo in range(self.ret_num_mo()):
            if eneocc:
                occ = self.ens[imo]
            else:
                occ = self.occs[imo]
            if abs(occ) < occmin: continue
            if abs(occ) > occmax: continue

            outstr += ' Sym= %s\n'%sym
            outstr += ' Ene= %f\n'%self.ens[imo]
            outstr += ' Spin= Alpha\n'
            outstr += ' Occup= %f\n'%occ
            for ibf, coeff in enumerate(self.mo_mat[:,imo]):
                outstr += '%10i  % 10E\n'%(ibf+1, coeff)

        return outstr

    def read(self, lvprt=1, spin=0):
        """
        Read in MO coefficients from a molden File.
        spin: 0 - read all coefficients, 1 - only alpha, -1 - only beta
        """
        MO = False
        GTO = False
        ATOMS = False
        coor_unit = 1.
        mo_vecs = []
        mo_ind = 0
        self.syms = [] # list with the orbital descriptions. they are entered after Sym in the molden file.
        self.occs = [] # occupations
        self.ens  = [] # orbital energies (or whatever is written in that field)

        num_bas={'s':1,'p':3,'sp':4,'d':6,'f':10,'g':15}

        orient={'s':['1'],
                'p':['x','y','z'],
               'sp':['1','x','y','z'],
                'd':['x2','y2','z2','xy','xz','yz'],
                'f':['xxx', 'yyy', 'zzz', 'xyy', 'xxy', 'xxz', 'xzz', 'yzz', 'yyz', 'xyz'],
                'g':15*['?']}

        num_orb=0
        curr_at=-1
        tmp_data = [None, None, None]
        spin_flag = True

        self.header = ''

        fileh = open(self.file, 'r')

        fstr = fileh.read()
        if ('[D5' in fstr) or ('[5D' in fstr):
            num_bas['d']=5
            orient['d']=['D0', 'D+1', 'D-1', 'D+2', 'D-2']
        if ('F7]' in fstr) or ('7F]' in fstr):
            num_bas['f']=7
            orient['f']=['F0', 'F+1', 'F-1', 'F+2', 'F-2', 'F+3', 'F-3']
        if ('9G]' in fstr) or ('9G]' in fstr):
            num_bas['g']=9
            orient['g']=9*['?']

        fileh.seek(0) # rewind the file

        while True:
            try:
                line = next(fileh)
            except StopIteration:
                print("Finished parsing %s"%fileh.name)
                break

            words = line.replace('=',' ').split()

            if 'molden format' in line.lower():
                if not '[Molden Format]' in line:
                    print(" WARNING: the header may not be understood by Jmol:")
                    print(line, end=' ')
                    print(" This has to be changed to:")
                    print(" [Molden Format]")

                    line = '[Molden Format]'

            # what section are we in
            if '[' in line:
                MO = False
                GTO = False
                ATOMS = False

            if '[MO]' in line:
                if lvprt >= 2: print("Found [MO] tag")
                MO = True
                GTO = False
            # extract the information in that section
            elif MO:
                if '=' in line:
                    if 'ene' in line.lower():
                        tmp_data[0] = float(words[-1])
                    elif 'sym' in line.lower():
                        tmp_data[1] = words[-1]
                    elif 'occ' in line.lower():
                        tmp_data[2] = float(words[-1])
                    elif 'spin' in line.lower():
                        tmp_spin = words[-1].lower()
                        if spin == -1 and tmp_spin == 'alpha':
                            spin_flag = False
                        elif spin == 1 and tmp_spin == 'beta':
                            spin_flag = False
                        else:
                            spin_flag = True
                if not '=' in line:
                    try:
                        tmp_vec = [float(words[1])]
                    except:
                        if words==[]: break # stop parsing the file if an empty line is found

                        print(" ERROR in lib_mo, parsing the following line:")
                        print(line)
                        raise

                    for ibas in range(num_orb-1):
                        line = next(fileh)
                        words = line.split()
                        tmp_vec.append(float(words[1]))

                    if spin_flag:
                        mo_ind += 1
                        mo_vecs.append(tmp_vec)
                        self.ens.append(tmp_data[0])
                        self.syms.append(tmp_data[1])
                        self.occs.append(tmp_data[2])

            elif ('[GTO]' in line):
                GTO = True
                # extract the information in that section
            elif GTO:
                if len(words)==0: # empty line: atom is finished
                    curr_at = -1
                elif curr_at==-1:
                    curr_at = int(words[0])
                    self.num_at = max(curr_at, self.num_at)
                elif (len(words) >= 2) and (words[0].lower() in num_bas):
                  orbsymb = words[0].lower()

                  for i in range(num_bas[orbsymb]):
                    self.basis_fcts.append(basis_fct(curr_at, orbsymb, orient[orbsymb][i]))
                    label = self.basis_fcts[-1].label()
                    if not label in self.bf_labels:
                        self.bf_labels.append(label)

                    num_orb+=1

            elif '[atoms]' in line.lower():
                ATOMS = True
                if 'au' in line.lower():
                    coor_unit = units.length['A']
            elif ATOMS:
                words = line.split()
                self.at_dicts.append({'Z':int(words[2]), 'x':float(words[3])*coor_unit, 'y':float(words[4])*coor_unit, 'z':float(words[5])*coor_unit})

            if not MO:
                self.header += line

        fileh.close()

### file parsing finished ###

        if lvprt >= 1 or len(mo_vecs[0])!=num_orb:
            print('\nMO file %s parsed.'%self.file)
            print('Number of atoms: %i'%self.num_at)
            print('Number of MOs read in: %i'%len(mo_vecs))
            print('Dimension: %i,%i,...,%i'%(len(mo_vecs[0]),len(mo_vecs[1]),len(mo_vecs[-1])))
            print('Number of basis functions parsed: ', num_orb)

        if len(mo_vecs[0])!=num_orb:
            raise error_handler.MsgError('Inconsistent number of basis functions!')

        if len(mo_vecs) > 1.8 * num_orb:
            print("""\n   WARNING: There are twice as many MOs as basis functions!
   If this is an unrestricted calculation use analyze_tden_unr.py\n""")

        if len(mo_vecs[-1]) == 0:
            lv = [len(mo_vec) for mo_vec in mo_vecs]
            imax = lv.index(0)
            print('*** WARNING: MO file contains MO vectors of zero length!')
            print('Using only the first %i entries'%imax)
            mo_vecs = mo_vecs[:imax]

        try:
           self.mo_mat = numpy.array(mo_vecs).transpose()
        except ValueError:
           print("\n *** Unable to construct MO matrix! ***")
           print("Is there a mismatch between spherical/cartesian functions?\n ---")
           raise

class MO_set_tddftb(MO_set):
    def read(self, lvprt=1, spin=0):
        """
        Read in MO coefficients from eigenvec.out file.
        Author: Ljiljana Stojanovic
        """
        if not spin==0:
            raise error_handler.MsgError("Only spin=0 implemented")

        MO = False
        mo_vecs = []
        mo_ind = 0
        mo_mat = []
        mo_array = []
        self.syms = [] # list with the orbital descriptions. they are entered after Sym in the molden file.
        self.occs = [] # occupations
        self.ens  = [] # orbital energies (or whatever is written in that field)

       # self.header = ''

        filemo = open('eigenvec.out', 'r')
        MO = False

        for line in filemo:
            words = line.split()
            if 'Eigenvector' in line:
                MO = True
                mo_vecs.append([])
            elif MO:
                 mo_ind += 1
                 if len(words) != 0:
                    if len(words) == 3:
                       mo_vecs[-1].append(float(words[1]))
                    elif len(words) == 5:
                       mo_vecs[-1].append(float(words[3]))
        filemo.close()

        try:
           mo_array = numpy.array(mo_vecs)
           self.mo_mat = numpy.array(mo_array).transpose()
        except ValueError:
           print("\n *** Unable to construct MO matrix! ***")
           print("Is there a mismatch between spherical/cartesian functions?\n ---")
           raise

        Eig = False
        Occs = False

        filemosec = open('detailed.out', 'r')
        for line in filemosec:
            if 'Eigenvalues /eV' in line:
               Eig = True
            elif 'Fillings' in line:
               Occs = True
            elif line.strip() == '':
               Eig = False
               Occs = False
            elif Eig:
                self.ens.append(float(line))
                self.syms.append(str('a'))
            elif Occs:
                self.occs.append(float(line))
        filemosec.close()

        if len(self.occs) == 0:
            self.ens = []
            self.syms = []
            print(" WARNING: Parsing of detailed.out failed - using band.out instead")
            ft = open('band.out', 'r')
            for line in ft.readlines()[1:]:
                words = line.split()
                if len(words) == 0:
                    break
                self.ens.append(float(words[1]))
                self.occs.append(float(words[2]))
            ft.close()

        at_symb = []

        filegeom = open('geom.xyz','r')
        nline = 0
        for line in filegeom:
            nline += 1
            words = line.split()
            if nline > 2:
               self.at_dicts.append({'Z':'', 'x':words[1], 'y':words[2], 'z':words[3]})
               at_symb.append(str(words[0]))
               self.num_at += 1
        filegeom.close()

        num_bas={'s':1,'sp':3,'spd':5}
        ang_mom={0:'s',1:'p',2:'d'}

        orient={'s':['1'],
                'sp':['1','2','3'],
                'spd':['1','2','3','4','5']}

        num_orb = 0
        ang_momentum = 0
        curr_at = 0
        for i in range(0,self.num_at):
            filewfc = open('wfc.3ob-3-1.hsd','r')
            atom_name = at_symb[i] + " {"
            for line in filewfc:
                atom = False
                words = line.split()
                if atom_name in line:
                 curr_at += 1
                 for line in filewfc:
                       words = line.split()
                       if (len(words) == 0): break
                       if 'AngularMomentum' in line:
                               ang_momentum = int(words[2])
                               if ang_momentum == 0:
                                    orbsymb = 's'
                                    orb_deg = 1
                                    num_orb = num_orb + orb_deg
                               elif ang_momentum == 1:
                                    orbsymb = 'sp'
                                    orb_deg = 3
                                    num_orb = num_orb + orb_deg
                               elif ang_momentum == 2:
                                    orbsymb = 'spd'
                                    orb_deg = 5
                                    num_orb = num_orb + orb_deg
                               for i in range(num_bas[orbsymb]):
                                   self.basis_fcts.append(basis_fct(curr_at, orbsymb, orient[orbsymb][i]))
                                   label = self.basis_fcts[-1].label()
                                   if not label in self.bf_labels:
                                      self.bf_labels.append(label)
            filewfc.close()

### file parsing finished ###

        if lvprt >= 1 or len(mo_vecs[0])!=num_orb:
            print('\nMO file %s parsed.'%self.file)
            print('Number of atoms: %i'%self.num_at)
            print('Number of MOs read in: %i'%len(mo_vecs))
            print('Dimension: %i,%i,...,%i'%(len(mo_vecs[0]),len(mo_vecs[1]),len(mo_vecs[-1])))
            print('Number of basis functions parsed: ', num_orb)

        if len(mo_vecs[0])!=num_orb:
            raise error_handler.MsgError('Inconsistent number of basis functions!')

class MO_set_adf(MO_set):
    """
    MO_set class for ADF.
    Note that ADF uses STOs!
    """
    def read(self, lvprt=1):
        """
        Extract the MOs from the TAPE21 file.
        Initial code written by S. Mai.
        """
        try:
            from scm.plams import KFFile
        except ImportError:
            from kf import kffile as KFFile

        f = KFFile(self.file)

        try:
            self.num_at = int(f.read('Geometry','nr of atoms'))
        except:
            print(("\n  ERROR: reading TAPE21 file (%s)!\n"%self.file))
            raise

        natomtype=int(f.read('Geometry','nr of atomtypes'))
        atomorder=f.read('Geometry','atom order index')
        atomindex=f.read('Geometry','fragment and atomtype index')
        atcharge=f.read('Geometry', 'atomtype total charge')
        nbptr=f.read('Basis','nbptr')
        naos2=int(f.read('Basis','naos'))

        # get number of basis functions per atomtype
        nbasis={}
        for iaty in range(natomtype):
            nbasis[iaty]=nbptr[iaty+1]-nbptr[iaty]

        if lvprt >= 2:
            print('Number of basis functions per atomtype:')
            print(nbasis)

        # get which atom is of which atomtype
        atom_type={}
        for iatom in range(self.num_at):
            index=atomorder[iatom]-1
            atom_type[iatom]=atomindex[index+self.num_at]-1

        if lvprt >= 2:
            print('Mapping of atoms on atomtypes:')
            print(atom_type)

        # get number of basis functions
        naos=0
        for i in atom_type:
            naos+=nbasis[atom_type[i]]
        if lvprt >= 2:
            print('Total number of basis functions:', naos)
        assert naos==naos2, "wrong number of orbitals"

        # map basis functions to atoms
        for iatom in range(self.num_at):
          index=atomorder[self.num_at+iatom]-1
          nb=nbasis[atom_type[index]]
          for iao in range(nb):
            self.basis_fcts.append(basis_fct(index+1))

        if lvprt >= 3:
            print('Basis functions:')
            for bf in self.basis_fcts:
                print(bf)

        # get MO coefficients
        NAO=int(f.read('Basis','naos'))
        npart=f.read('A','npart')
        NMO_A=int(f.read('A','nmo_A'))
        mocoef_A=f.read('A','Eigen-Bas_A')

        for i in range(len(npart)):
          npart[i]-=1
        self.mo_mat = numpy.zeros([NAO, NMO_A])
        iao=0
        imo=0
        for i,el in enumerate(mocoef_A):
          iao1=npart[iao]
          self.mo_mat[iao1][imo]=el
          iao+=1
          if iao>=NAO:
            iao=0
            imo+=1

        nelec = int(f.read('General','electrons'))
        assert nelec%2==0, "Odd number of electrons not supported"
        nocc = nelec // 2
        nvirt = NMO_A - nocc
        self.occs = nocc * [2.] + nvirt * [0.]

        # read coordinates and elements
        atom_Z=[int(atcharge[atom_type[i]]) for i in range(self.num_at)]
        xyz = numpy.array(f.read('Geometry', 'xyz')).reshape(self.num_at, 3) * units.length['A']

        for iat in range(self.num_at):
            self.at_dicts.append({'Z':atom_Z[iat], 'x':xyz[atomorder[iat]-1, 0], 'y':xyz[atomorder[iat]-1, 1], 'z':xyz[atomorder[iat]-1, 2]})

## get AO overlap matrix
## NOTE: Smat is only read from TAPE15 !
## hence, if necessary, say "SAVE TAPE21 TAPE15" in input
#if f2:
#  NAO = f2.read('Basis','naos')
#  Smat = f2.read('Matrices','Smat')
#
#  ao_ovl=[ [ 0. for i in range(NAO) ] for j in range(NAO) ]
#  x=0
#  y=0
#  for el in Smat:
#    ao_ovl[x][y]=el
#    ao_ovl[y][x]=el
#    x+=1
#    if x>y:
#      x=0
#      y+=1

class MO_set_onetep(MO_set):
    """
    Handling of MOs for ONETEP.
    ONETEP actually already provides the 1TDM already in the AO basis.
    That is why the workflow of the ONETEP interface is different to the others.
    """
    def read(self, lvprt=1):

        self.num_at = 0
        if lvprt >= 1:
            print("Reading basis function information ...")
        for line in open(self.file + '.jointngwf2atoms'):
            if line[0] == '#':
                continue
            at_ind = int(line.split()[-1])
            self.basis_fcts.append(basis_fct(at_ind))
            self.num_at = max(self.num_at, at_ind)

        # The MO coefficients are just a dummy unit matrix.
        self.mo_mat = numpy.identity(len(self.basis_fcts))

        # Information for valence set only
        self.basis_fcts_val = []
        for line in open(self.file + '.valngwf2atoms'):
            if line[0] == '#':
                continue
            at_ind = int(line.split()[-1])
            self.basis_fcts_val.append(basis_fct(at_ind))
        assert(len(self.basis_fcts_val) == self.ret_num_bas() // 2)

        if lvprt >= 1:
            print("Reading and processing overlap matrix ...")
        self.read_S()

    def read_S(self):
        """
        Read the overlap matrix.
        """
        # analysis of joint set
        # this uses the variables of the main class
        tmp = []
        Sname = self.file + '.jointoverlap.mat'
        for line in open(Sname, 'r'):
            tmp.append(line.split())
        self.S = numpy.array(tmp, float)
        (eval, evec) = numpy.linalg.eigh(self.S)
        # sqval = numpy.sqrt(eval)
        # self.S2 = numpy.dot(evec,
        #     numpy.dot(numpy.diag(sqval), evec))

        # analysis of valence set only
        # this creates its own temporary variables
        tmp = []
        Sname = self.file + '.valoverlap.mat'
        for line in open(Sname, 'r'):
            tmp.append(line.split())
        self.S_val = numpy.array(tmp, float)

    def lowdin_trans(self, D):
        """
        Nothing to do. D is already stored in the Lowdin orthogonalized AO basis.
        """
        raise error_handler.MsgError("Not implemented for ONETEP")

    def OmBas_Mulliken(self, D, formula):
        DS = numpy.dot(D, self.S)
        SD = numpy.dot(self.S_val, D)
        if   formula == 0:
            return DS * SD
        elif formula == 1:
            SDS = numpy.dot(SD, self.S)
            return 0.5 * (DS * SD + D * SDS)
        else:
            raise error_handler.ElseError(formula, 'Om formula')

    def comp_OmAt(self, OmBas):
        """
        Compute the Omega matrix wrt atoms.
        Differentiate between valence and joint basis set.
        """
        OmAt = numpy.zeros([self.num_at, self.num_at])
        for iat, ist, ien in self.bf_blocks(self.basis_fcts_val, self.ret_num_bas() // 2):
            for jat, jst, jen in self.bf_blocks():
                OmAt[iat, jat] = numpy.sum(OmBas[ist:ien, jst:jen])

        return OmAt

class basis_fct:
    """
    Container for basisfunction information.
    """
    def __init__(self, at_ind=-1, l='?', ml='?'):
        self.set(at_ind, l, ml)

    def set(self, at_ind=-1, l='?', ml='?'):
        self.at_ind = at_ind # atom where the function is located, starting at 1
        self.l = l   # s, p, d, f
        self.ml = ml # x, y, z

    def __str__(self):
        return "Basis function (%s, %s) at atom %i"%(self.l, self.ml, self.at_ind)

    def label(self):
        return "%s-%s"%(self.l, self.ml)

class jmol_MOs:
    """
    Class for producing input for the Jmol program that can be used to plot MOs.
    """
    def __init__(self, name, width=400):
        self.name = name
        self.imo = -1
        self.width = width

    def pre(self, ofile=None):
        self.jmfile = open("%s_jmol.spt"%self.name,'w')
        if not ofile==None:
            self.jmfile.write('load %s FILTER "nosort"\n'%ofile)
        self.jmfile.write ('mo titleformat ""\n\n')
        self.jmfile.write('background white\nmo fill\nmo cutoff 0.04\n')

        self.htmlfile = lib_file.htmlfile("%s.html"%self.name)
        self.htmlfile.pre(self.name)

    def add_mo(self, jmolstr, name, val):
        imfile = "%s_%.2f.png"%(name,abs(val))

        self.jmfile.write(jmolstr)
        self.jmfile.write('write image pngt "%s"\n'%imfile)

        if self.imo%2==0: self.htmlfile.write("<tr>\n")
        self.htmlfile.write("<td>")
        self.htmlfile.write("<img src=\"%s\" "%(imfile))
        self.htmlfile.write("border=\"1\" width=\"%i\">"%self.width)
        self.htmlfile.write("<br> %.3f"%val)
        self.htmlfile.write("</td>\n")
        if self.imo%2==1: self.htmlfile.write("</tr>\n")

        self.imo += 1

    def next_set(self,header):
        if self.imo%2==1: self.htmlfile.write("</tr>\n")
        if not self.imo==-1: self.htmlfile.write("</table>\n")

        self.imo = 0

        self.jmfile.write("\n")

        self.htmlfile.write("\n<h2>%s</h2>\n"%header)
        self.htmlfile.write("<table>\n")

    def post(self):
        self.jmfile.close()

        if self.imo%2==1: self.htmlfile.write("</tr>\n")

        self.htmlfile.write("</table>\n")
        self.htmlfile.post()

        print("\nJmol input file %s and %s written"%(self.jmfile.name, self.htmlfile.name))
        print("   Run as: jmol -n %s"%self.jmfile.name)
