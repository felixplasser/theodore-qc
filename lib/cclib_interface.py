"""
This is an interface to the external cclib Library.
http://cclib.github.io
Download and install cclib if you want to use the functions.
"""

import struct
import numpy
import file_parser, lib_mo, error_handler, units, lib_struc
try:
    import openbabel
except ImportError:
    print(" *** Warning: python-openbabel not found! ***")
    print(" Using emulation program with limited capabilities ...")
    import OB_repl as openbabel

class file_parser_cclib(file_parser.file_parser_base):
    def __init__(self, *args, **kwargs):
        file_parser.file_parser_base.__init__(self, *args, **kwargs)

        self.data, self.prog = self.get_data(self.ioptions['rtype'], self.ioptions['rfile'])

    def get_data(self, rtype, rfile):
        """
        Call cclib to read in the data.
        """
        try:
            import cclib.parser
            print("Using cclib version %s"%cclib.__version__)
        except ImportError:
            print("\n ERROR: Did not find the external cclib package!")
            print("Please, download and install cclib to use this functionality.\n")
            raise

        cfile = cclib.parser.ccopen(rfile)
        if cfile == None:
            raise error_handler.MsgError('File %s cannot be parsed by cclib!'%rfile)
        else:
            print("Parsing %s ..."%str(cfile))

        prog = str(cfile).split()[0]

        return cfile.parse(), prog

    def read(self, mos, rect_dens=True):
        """
        Convert the information to data analyzed by TheoDORE.
        """
        state_list = []

        try:
            core_orbs = sum(self.data.coreelectrons) / 2
            if core_orbs > 0:
                print("\n WARNING: %i core orbitals detected"%core_orbs)
                print(" Please check the results carefully!\n")
                rect_dens = False
        except AttributeError:
            pass

        for ist, exc_en in enumerate(self.data.etenergies):
            state_list.append({})
            state = state_list[-1]

            state['exc_en'] = exc_en / units.energy['rcm'] * units.energy['eV']
            try:
                state['osc_str'] = self.data.etoscs[ist]
            except AttributeError:
                state['osc_str'] = -1.
            state['irrep'] = self.data.etsyms[ist].replace('let-', '')

        if self.prog == 'ADF':
            self.tden_adf(state_list, mos, rect_dens)
        elif self.prog == 'ORCA':
            if self.ioptions['read_binary']:
                self.tden_orca(state_list, mos, rect_dens)
            else:
                print 'Reading ORCA input from stdout (read_binary=False)'
                self.tden_cclib(state_list, mos, rect_dens)
        else:
            self.tden_cclib(state_list, mos, rect_dens)

        return state_list

    def tden_cclib(self, state_list, mos, rect_dens):
        for ist, state in enumerate(state_list):
            state['state_ind'] = ist + 1
            state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
            print("\n" + state['name'])

            state['tden'] = self.init_den(mos, rect=rect_dens)
            for conf in self.data.etsecs[ist]:
                [(iocc, spocc), (ivirt, spvirt), coeff] = conf
                assert(iocc < ivirt) # de-excitations cannot be handled at this point
                state['tden'][iocc, ivirt] = coeff

                if coeff*coeff > 0.05:
                    print(" (%i->%i), coeff=% .4f"%(iocc+1, ivirt+1, coeff))
                elif numpy.isnan(coeff):
                    print(" ERROR: (%i->%i), coeff=% .4f"%(iocc+1, ivirt+1, coeff))
                    if self.prog == 'ORCA':
                        print("For ORCA/RPA jobs, please set read_binary=True\n")
                    raise error_handler.MsgError("Not-a-number encountered")

    def tden_adf(self, state_list, mos, rect_dens):
        print("\n   WARNING: using deprecated old ADF interface!\n")
        import kf

        print("Opening file TAPE21 ...")
        rfile = kf.kffile('TAPE21')

        try:
            nmo = int(rfile.read('A','nmo_A'))
            nelec = int(rfile.read('General','electrons'))
        except:
            print("\n  ERROR: reading file TAPE21!\n")
            raise
        assert nelec%2==0, "Odd number of electrons"

        group = rfile.read('Geometry','grouplabel')[0]
        if not group in ['NOSYM', 'C1', 'c1']:
            print("grouplabel: " + group)
            raise error_handler.MsgError('No support for symmetry')

        nocc = nelec / 2
        assert nocc == mos.ret_ihomo() + 1
        nvirt = nmo - nocc
        assert nvirt == mos.ret_num_mo() - nocc

        try:
            nsing = len(rfile.read('All excitations','All Sing-Sing excitations'))
        except TypeError:
            nsing = 0
        try:
            ntrip = len(rfile.read('All excitations','All Sing-Trip excitations'))
        except TypeError:
            ntrip = 0

        assert nsing+ntrip == len(state_list)

        istate = 0
        for ising in range(nsing):
            state = state_list[istate]

            state['state_ind'] = ising + 1
            state['mult'] = 1
            state['name'] = '%i(%i)%s'%(state['state_ind'], state['mult'] ,state['irrep'])

            state['tden'] = self.init_den(mos, rect=rect_dens)
            eigen = rfile.read('Excitations SS A','eigenvector %s'%(istate+1))
            state['tden'][:,nocc:nmo] = eigen.reshape(nocc, nvirt)
            istate += 1

            # print-out
            print(state['name'])
            tden = state['tden']
            for i in range(len(tden)):
                for j in range(len(tden[0])):
                    val = tden[i, j]
                    if val*val > 0.1:
                        print("(%i -> %i) % .4f"%(i+1,j+1,val))

        for itrip in range(ntrip):
            state = state_list[istate]

            state['state_ind'] = itrip + 1
            state['mult'] = 3
            state['name'] = '%i(%i)%s'%(state['state_ind'], state['mult'] ,state['irrep'])

            state['tden'] = self.init_den(mos, rect=rect_dens)
            eigen = rfile.read('Excitations ST A','eigenvector %s'%(istate+1 - nsing))
            state['tden'][:,nocc:nmo] = eigen.reshape(nocc, nvirt)
            istate += 1

            # print-out
            print(state['name'])
            tden = state['tden']
            for i in range(len(tden)):
                for j in range(len(tden[0])):
                    val = tden[i, j]
                    if val*val > 0.1:
                        print("(%i -> %i) % .4f"%(i+1,j+1,val))

        rfile.close()

    def tden_orca(self, state_list, mos, rect_dens, filen='orca.cis'):
        """
        Read binary CI vector file from ORCA.
        Authors: S. Mai, F. Plasser
        """
        print "Reading CI vectors from binary ORCA file %s"%filen

        CCfile=open(filen,'rb')
        # header
        # consists of 9 4-byte integers, the first 5 of which give useful info
        nvec  =struct.unpack('i', CCfile.read(4))[0]
        print "Number of vectors:", nvec
        header=[ struct.unpack('i', CCfile.read(4))[0] for i in range(8) ]
        #print header

        # header array contains:
        # [0] index of first alpha occ,  is equal to number of frozen alphas
        # [1] index of last  alpha occ
        # [2] index of first alpha virt
        # [3] index of last  alpha virt, header[3]+1 is equal to number of bfs
        # [4] index of first beta  occ,  for restricted equal to -1
        # [5] index of last  beta  occ,  for restricted equal to -1
        # [6] index of first beta  virt, for restricted equal to -1
        # [7] index of last  beta  virt, for restricted equal to -1

        if any( [ header[i]!=-1 for i in range(4,8) ] ):
            raise error_handler.MsgError("No support for nrestricted MOs")

        nfrzc = header[0]
        nocc = header[1] + 1
        nact = nocc - nfrzc
        nmo  = header[3] + 1
        nvir = nmo - header[2]
        lenci = nact*nvir
        print '  nmo: %i , nocc: %i , nact: %i , nvir: %i'%(nmo,nocc,nact,nvir)

        # loop over states
        # for non-TDA order is: X+Y of 1, X-Y of 1, X+Y of 2, X-Y of 2, ...
        prevroot = -1
        istate = 0
        for ivec in range(nvec):
            # header of each vector
            # contains 6 4-byte ints, then 1 8-byte double, then 8 byte unknown
            nele,d1,mult,d2,iroot,d3 = struct.unpack('iiiiii', CCfile.read(24))
            ene,d3 = struct.unpack('dd', CCfile.read(16))
            print '  nele: %i , mult: %i , iroot: %i'%(nele,mult,iroot)
            #print d1,d2,d3
            #print '  Ene: %.4f, Osc: %.4f'%(ene*units.energy['eV'], osc)
            # -> energy does not look consistent

            # then comes nact * nvirt 8-byte doubles with the coefficients
            coeff = struct.unpack(lenci*'d', CCfile.read(lenci*8))

            if prevroot!=iroot:
                state = state_list[istate]
                state['state_ind'] = iroot + 1
                state['mult'] = mult
                state['irrep'] = 'A'
                state['name'] = '%i(%i)%s'%(state['state_ind'], state['mult'] ,state['irrep'])
                state['tden'] = self.init_den(mos, rect=rect_dens)
                #assert abs(state['exc_en']-ene*units.energy['eV']) < 1.e-4, 'Incorrect energies'
                # -> does not work for all states

                state['tden'][nfrzc:nocc, nocc:nmo] = numpy.reshape(coeff, [nact,nvir])
                istate += 1
            else:
            # in this case, we have a non-TDA state!
            # and we need to compute (prevvector+currentvector)/2 = X vector
                print 'Constructing X-vector of RPA state'
                state['tden'][nfrzc:nocc, nocc:nmo] += numpy.reshape(coeff, [nact,nvir])
                state['tden'][nfrzc:nocc, nocc:nmo] *= .5

            prevroot=iroot

    def check(self, lvprt=1):
        """
        Check if the input file can be used by cclib.
        """
        errcode = 0

        if lvprt >= 1:
            print("\nChecking if the logfile %s can be parsed by cclib ..."%self.ioptions['rfile'])
            print("\nEssential attributes:")

        for attr in ['mocoeffs', 'atombasis', 'natom', 'homos', 'moenergies', 'etenergies', 'etsyms', 'etsecs']:
            chk = hasattr(self.data, attr)
            if not chk: errcode = 2

            if lvprt >= 1:
                print('%15s ... %s'%(attr, chk))

        if lvprt >= 1:
            print("\nOptional attributes:")

        for attr in ['etoscs', 'aooverlaps', 'mosyms']:
            chk = hasattr(self.data, attr)

            if lvprt >= 1:
                print('%15s ... %s'%(attr, chk))

        if lvprt >= 1:
            print("\nAttributes for structure parsing and creation of Molden file:")

        for attr in ['gbasis', 'natom', 'atomcoords', 'atomnos']:
            chk = hasattr(self.data, attr)
            if not chk: errcode = max(1, errcode)

            if lvprt >= 1:
                print('%15s ... %s'%(attr, chk))

        return errcode

    def read_mos(self):
        """
        Return an MO_set object in TheoDORE format.
        """
        if self.prog == 'ADF':
            #mos = MO_set_adf(file = None)
            mos = MO_set_cclib(file = None)
        else:
            mos = MO_set_cclib(file = None)

        mos.read(self.data, lvprt=2)

        return mos

    def ret_struc(self,lvprt=1):
        struc = structure_cclib()
        struc.read_cclib(self.data)

        return struc

class MO_set_cclib(lib_mo.MO_set_molden):
    def read(self, data, lvprt=1):
        """
        Read cclib data and convert to TheoDORE format.
        """
        self.mo_mat = data.mocoeffs[0].transpose()
        if lvprt >= 1:
            print("MO-matrix read from cclib, dimension: %i x %i"%(len(self.mo_mat), len(self.mo_mat[0])))

        self.num_at = data.natom

        self.ens = data.moenergies[0]
        self.ihomo = data.homos[0]
        self.occs = (self.ihomo + 1) * [1.] + (len(self.ens) - self.ihomo - 1) * [0.]

        self.basis_fcts = [lib_mo.basis_fct(-1) for ibas in range(self.ret_num_bas())]
        maxbas = -1
        for iat, baslist in enumerate(data.atombasis):
            maxbas = max([maxbas] + baslist)
            for ibas in baslist:
                self.basis_fcts[ibas].set(iat + 1)

        assert(self.ret_num_bas() == maxbas+1)

        try:
            self.syms = data.mosyms
        except AttributeError:
            self.syms = ['X'] * len(self.ens)

        try: self.S = data.aooverlaps
        except AttributeError: pass

        try: self.atomnos = data.atomnos
        except AttributeError: pass

        try: self.coords = data.atomcoords[-1]
        except AttributeError: pass

        try: self.gbasis = data.gbasis
        except AttributeError: pass

    def ret_ihomo(self):
        return self.ihomo

    def write_molden_file(self, fname='out.mld', cfmt='% 10E', occmin=-1):
        """
        Write a file in Molden format.
        """
        self.header = "[Molden Format]\n"

        self.header+= "[Atoms] Angs\n"
        for iat in xrange(self.num_at):
            atno = self.atomnos[iat]
            (x, y, z) = self.coords[iat]
            self.header+= "%3s %5i %3i %12.6f %12.6f %12.6f\n"%\
                (lib_struc.Z_symbol_dict[atno], iat+1, atno, x, y, z)

        self.header+= "[GTO]\n"
        for iat, at in enumerate(self.gbasis):
            self.header+= "%3i 0\n"%(iat+1)
            for bas in at:
                typ = bas[0]
                exps = bas[1]

                self.header+= "%3s %5i 1.00\n"%(typ, len(exps))

                for exp in exps:
                    self.header+= "%18.10E %18.10E\n"%exp
            self.header+= "\n"

        self.export_AO(self.ens, self.occs, self.mo_mat.transpose(), fname, cfmt, occmin)

class MO_set_adf(lib_mo.MO_set_molden):
    def read(self, data, lvprt=1):
        """
        Read MOs from TAPE21 file.
        """
        print("\n   Using new ADF interface\n")
        import kf

        print("Opening file TAPE21 ...")
        rfile = kf.kffile('TAPE21')

        NAO = rfile.read('Basis','naos')
        NMO = rfile.read('A','nmo_A')

        self.mo_mat = rfile.read('A','Eigen-Bas_A').reshape(NAO, NMO)
        if lvprt >= 1:
            print("MO-matrix read from TAPE21, dimension: " + str(self.mo_mat.shape))

        self.num_at = data.natom

        self.ens = rfile.read('A','eps_A')
        self.ihomo = data.homos[0]
        self.occs = rfile.read('A','froc_A')

        npart = rfile.read("A","npart")

        print npart
        print data.atombasis

        # assign basis functions to atoms
        # This can be found at
        # BAS: List of all Elementary Cartesian Basis Functions

        self.basis_fcts = [lib_mo.basis_fct(-1) for ibas in range(self.ret_num_bas())]
        maxbas = -1
        for iat, baslist in enumerate(data.atombasis):
            maxbas = max([maxbas] + baslist)
            print baslist
            for ibas in baslist:
                self.basis_fcts[ibas].set(iat + 1)

        assert(self.ret_num_bas() == maxbas+1)

        rfile.close()

        raise error_handler.NIError

class structure_cclib(lib_struc.structure):
    def read_cclib(self, data, ind=-1, lvprt=1):
        if lvprt>=1:
            print("Reading cclib structure with %i atoms."%data.natom)

        self.mol = openbabel.OBMol()

        for iat in xrange(data.natom):
            obatom = openbabel.OBAtom()
            obatom.SetAtomicNum(int(data.atomnos[iat]))
            coords = data.atomcoords[ind][iat]
            obatom.SetVector(*coords)

            self.mol.AddAtom(obatom)
