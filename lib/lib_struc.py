"""
Tools for molecular structure analysis and manipulation.
This is a wrapper to python-openbabel.
"""

import os, shutil, locale
import numpy
try:
    import openbabel
except ImportError:
    print " *** Warning: python-openbabel not found! ***"
    print " Using emulation program with limited capabilities ..."
    import OB_repl as openbabel
import units

Z_symbol_dict = {1:'H',5:'B',6:'C',7:'N',8:'O',9:'F',15:'P',16:'S',17:'Cl',26:'Fe',34:'Se',35:'Br',50:'Sn',53:'I',77:'Ir'}
symbol_Z_dict = {}
for key,val in Z_symbol_dict.iteritems():
    symbol_Z_dict[val] = key
Z_symbol_dict[99] = 'X'
veloc_types = ['vtxyz','vnx'] # these are defined below

class structure:
    """
    Class to manipulate a structure.
    """
    def __init__(self, name=''):
        self.name = name
        self.new_types = ['txyz2','col','colr','nx'] # these are defined here

    def read_file(self, file_path, file_type='tmol'):
        """
        Read in the structure from a file.
        """
        self.file_path = file_path
        self.file_type = file_type
        self.mol = openbabel.OBMol()

        if self.file_type in self.new_types:
            self.read_new_type()
        else:
            obconversion = openbabel.OBConversion()
            obconversion.SetInFormat(file_type)
            obconversion.ReadFile(self.mol, file_path)

    def read_new_type(self):
        """
        Routine for reading a file type not contained in open babel.
        """
        infile = open(self.file_path,'r')
        line = infile.readline()

        if self.file_type == 'txyz2': # Tinker format preserving all information
           self.tinker_symbs = []
           self.tinker_extra = []
           num_at_chk = int(line.split()[0])
           line = infile.readline()
           while(line!=''):
              words = line.split()
              if len(words) == 0:
                  break
              obatom = openbabel.OBAtom()
              self.tinker_symbs.append(words[1])
              lett1 = words[1][0]
              if words[1] in symbol_Z_dict:
                 at_num = symbol_Z_dict[words[1]]
              elif words[1][0] in symbol_Z_dict:
                 at_num = symbol_Z_dict[words[1][0]]
              else:
                 at_num = 99
              obatom.SetAtomicNum(at_num)

              coords = [float(word) for word in words[2:5]]
              obatom.SetVector(*coords)
              self.tinker_extra.append(words[5:])
              # maybe one could use SetSymbol here

              self.mol.AddAtom(obatom)
              line = infile.readline()

           #print self.mol.NumAtoms()
           assert(num_at_chk==self.mol.NumAtoms())
        elif self.file_type == 'col' or self.file_type == 'nx':
           while(line!=''):
              words = line.split()
              obatom = openbabel.OBAtom()
              obatom.SetAtomicNum(int(float(words[1])))
              coords = [float(word)*units.length['A'] for word in words[2:5]]
              obatom.SetVector(*coords)

              self.mol.AddAtom(obatom)
              line = infile.readline()
        else:
           print 'type %s not supported for input'%self.file_type
           exit(1)

        infile.close()

    def get_mol(self, mol, file_path, file_type='xyz'):
        """
        Read in an openbabel mol that is passed from a different routine.
        Can be used for accessing multiple structure xyz files.
        """
        self.file_path = file_path
        self.file_type = file_type
        self.mol = mol

    def read_file_vector(self, def_file_path, file_type, vector):
        """
        Initialise the structure with a default .mol from a file and a vector.
        """
        self.read_file(def_file_path, file_type) # like this because copying objects doesn't work

        for i in xrange(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            atom.SetVector(vector[3*i], vector[3*i+1], vector[3*i+2])
        #print 'read_file_vector done'

    def read_file_3xN_matrix(self, def_file_path, file_type, coor_mat):
        """
        Initialise the structure with a default .mol and a 3xN matrix.
        """
        self.read_file(def_file_path, file_type)

        #print self.mol.NumAtoms(), file_type

        #for i, atom in enumerate(openbabel.OBMolAtomIter(self.mol)): ## this only works with the new version

        self.read_3xN_matrix(coor_mat)

    def read_3xN_matrix(self, coor_mat, at_list=None):
        """
        Read in a 3xN matrix.
        If an optional <at_list> is specified, it is assumed that the coordinates
          specified correspond to these atoms and only these atoms are changed.
        """
        if at_list == None:
            at_list = [i+1 for i in range(self.mol.NumAtoms())]

        for imat,iat in enumerate(at_list):
            atom = self.mol.GetAtom(iat)
            atom.SetVector(coor_mat[imat][0], coor_mat[imat][1], coor_mat[imat][2])

    def read_at_dicts(self, at_dicts):
        """
        Create a mol from information specified in a list of dictionaries.
        [{'Z':, 'x':, 'y':, 'z':}, ...]
        """
        self.mol = openbabel.OBMol()
        
        for iat in xrange(len(at_dicts)):
            obatom = openbabel.OBAtom()
            obatom.SetAtomicNum(at_dicts[iat]['Z'])
            coords = (at_dicts[iat]['x'], at_dicts[iat]['y'], at_dicts[iat]['z'])
            obatom.SetVector(*coords)
            
            self.mol.AddAtom(obatom)
            
    def ret_vector(self):
        " All the coordinates in one vector "
        vec_list = []
        for i in xrange(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            vec_list += [atom.x(), atom.y(), atom.z()]

        return numpy.array(vec_list)

    def ret_3xN_matrix(self, at_list=None):
        """
        Return coordinates in a 3 x N matrix.
        If <at_list> is specified only the atoms with those indices are considered.
        """
        mat_list = []

        if at_list == None:
            at_list = [i+1 for i in range(self.mol.NumAtoms())]

        for i in at_list:
            atom = self.mol.GetAtom(i)
            mat_list += [[atom.x(), atom.y(), atom.z()]]

        return numpy.array(mat_list)


    def ret_moved_structure(self, add_vec, name=''):
        """
        Move the structure by <add_vec>.
        """
        if name == '': name = self.name

        coor_mat = self.ret_3xN_matrix()
        coor_mat += add_vec

        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, coor_mat)
        return ret_struc

    def ret_bond_length(self, i, j):
        """
        Return the distance between atoms indexed i and j.
        """

        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)

        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])

        return numpy.dot(pos_i - pos_j, pos_i - pos_j)**.5

    def ret_distance_matrix(self):
        """
        Return a matrix containing all the distances between atoms.
        """
        num_at = self.mol.NumAtoms()

        ret_mat = numpy.zeros([num_at, num_at])

        for iat in xrange(num_at):
            for jat in xrange(iat+1, num_at):
                bij = self.ret_bond_length(iat+1, jat+1)
                ret_mat[iat, jat] = bij
                ret_mat[jat, iat] = bij

        return ret_mat

    def ret_bend(self, i, j, k):
        """
        Return the bending angle between atoms indexed i, j, k.
        """

        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)
        OBAtom_k = self.mol.GetAtom(k)

        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])
        pos_k = numpy.array([OBAtom_k.x(), OBAtom_k.y(), OBAtom_k.z()])

        vec1 = pos_i - pos_j
        vec2 = pos_k - pos_j

        len_1 = numpy.sqrt(numpy.dot(vec1,vec1))
        len_2 = numpy.sqrt(numpy.dot(vec2,vec2))

        return numpy.arccos(numpy.dot(vec1, vec2) / (len_1*len_2)) / numpy.pi * 180

    def ret_tors(self, i, j, k, l):
        """
        Return the torsion angle between atoms indexed i, j, k, l.
        """
        # Dihedral angle computed according to (http://en.wikipedia.org/wiki/Dihedral_angle) to get the full 360 deg range.

        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)
        OBAtom_k = self.mol.GetAtom(k)
        OBAtom_l = self.mol.GetAtom(l)

        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])
        pos_k = numpy.array([OBAtom_k.x(), OBAtom_k.y(), OBAtom_k.z()])
        pos_l = numpy.array([OBAtom_l.x(), OBAtom_l.y(), OBAtom_l.z()])

        vec1 = pos_j - pos_i
        vec2 = pos_k - pos_j
        vec3 = pos_l - pos_k

        cross1 = numpy.cross(vec1, vec2)
        cross2 = numpy.cross(vec2, vec3)

        norm2 = numpy.sqrt(numpy.dot(vec2, vec2))
        dot1 = numpy.dot(vec1, cross2)
        dot2 = numpy.dot(cross1, cross2)

        return numpy.arctan2(norm2 * dot1, dot2) / numpy.pi * 180

    def ret_symbol(self, i):
        """
        Returns the symbol of atom i.
        """
        try:
            return Z_symbol_dict[self.mol.GetAtom(i).GetAtomicNum()]
        except KeyError:
            return 'X'

    def ret_num_at(self):
        """
        Returns the number of atoms in the file.
        """
        return self.mol.NumAtoms()

    def ret_mass_vector(self, power):
        """
        Returns a vector with the masses of the atoms (each 1 time) taken to the <power> power.
        """
        mass_list = []
        for i in xrange(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            mass_list += [atom.GetAtomicMass()**power]

        return numpy.array(mass_list, float)

    def ret_partition(self,cutBonds=[], lvprt=1):
        """
        Return a partition according to different non-bonded molecules.
        cutBonds is a list of tuples for bonds to cut. The order does not matter
        e.g. cutBonds=[(3,4),(7,8)]
        """
        at_lists = []
        chk_list = []
        remaining_atoms = [self.mol.NumAtoms()-i for i in xrange(self.mol.NumAtoms())]

        while(len(remaining_atoms)>0): # do the loop as long as there are still atoms which are not in at_lists
            if len(chk_list) > 0:
                curr_at = chk_list.pop()
            else:
                at_lists.append([])
                curr_at = remaining_atoms.pop()
                at_lists[-1].append(curr_at)

            atom = self.mol.GetAtom(curr_at)
            for bonded in openbabel.OBAtomAtomIter(atom):
                bind = bonded.GetIdx()
                if bind in remaining_atoms:
                    if  (curr_at,bind) in cutBonds or (bind,curr_at) in cutBonds:
                        print 'cutting bond %i-%i'%(curr_at,bind)
                    else:
                        del remaining_atoms[remaining_atoms.index(bind)]
                        chk_list.append(bind)
                        at_lists[-1].append(bind)

        if lvprt >= 1:
            print "\n*** Fragment composition ***"
            for i, at_list in enumerate(at_lists):
                print "  Fragment %i: %s"%(i+1, self.ret_at_list_composition(at_list))

        return at_lists

    def ret_at_list_composition(self, at_list):
        """
        Return a string describing the atoms contained in at_list,
           e.g. C4H5N3
        """
        symb_list = []
        at_dict = {}

        for iat in at_list:
            symb = self.ret_symbol(iat)
            try:
                at_dict[symb] += 1
            except KeyError:
                at_dict[symb]  = 1
                symb_list.append(symb)

        ret_str = ''
        for symb in sorted(symb_list):
            numel = at_dict[symb]
            ret_str += '%s '%symb if numel==1 else '%s%i '%(symb, numel)

        return ret_str

    def make_coord_file(self, file_path, file_type='tmol'):
        """
        Write the structure file.
        """
        if file_type in self.new_types:
            self.make_coord_new(file_path, file_type)
        else:
            obconversion = openbabel.OBConversion()
            obconversion.SetOutFormat(file_type)
            obconversion.WriteFile(self.mol, file_path)

    def make_coord_new(self, file_path, file_type):
        """
        Routine for reading a file type not contained in open babel.
        """
        outfile = open(file_path,'w')
        num_at = self.mol.NumAtoms()

        if file_type == 'txyz2':
          outfile.write('%i from MSMT\n'%num_at)
          for ind in xrange(1, num_at+1):
            obatom = self.mol.GetAtom(ind)
            outstr  = ' %6i'%ind
            outstr += ' %4s'%self.tinker_symbs[ind-1]
            outstr += ' % 10.6f'%obatom.x()
            outstr += ' % 10.6f'%obatom.y()
            outstr += ' % 10.6f'%obatom.z()
            for extr in self.tinker_extra[ind-1]:
                outstr += ' %6s'%extr
            outfile.write(outstr+'\n')
        elif file_type == 'col':
          for ind in xrange(1, num_at+1):
            obatom  = self.mol.GetAtom(ind)
            outstr  = '%2s'%Z_symbol_dict[obatom.GetAtomicNum()]
            outstr += ' %7.1f'%obatom.GetAtomicNum()
            outstr += '% 14.8f'%(obatom.x() / units.length['A'])
            outstr += '% 14.8f'%(obatom.y() / units.length['A'])
            outstr += '% 14.8f'%(obatom.z() / units.length['A'])
            #outstr += '   %12.8f'%obatom.GetAtomicMass()
            # take the mass of the most abundant isotope rather than the average
            outstr += '% 14.8f'%obatom.GetExactMass()
            outfile.write(outstr+'\n')
        elif file_type == 'colr':
          print 'colr specified - atoms will be ordered by type'
          atnums=[]
          atstrs={}
          for ind in xrange(1, num_at+1):
            obatom  = self.mol.GetAtom(ind)

            outstr  = '%2s'%Z_symbol_dict[obatom.GetAtomicNum()]
            outstr += ' %7.1f'%obatom.GetAtomicNum()
            outstr += '% 14.8f'%(obatom.x() / units.length['A'])
            outstr += '% 14.8f'%(obatom.y() / units.length['A'])
            outstr += '% 14.8f'%(obatom.z() / units.length['A'])
            #outstr += '   %12.8f'%obatom.GetAtomicMass()
            # take the mass of the most abundant isotope rather than the average
            outstr += '% 14.8f\n'%obatom.GetExactMass()

            Z = obatom.GetAtomicNum()

            if not Z in atnums:
                atnums.append(Z)
                atstrs[Z] = []

            atstrs[Z].append(outstr)

          for Z in atnums:
              for ostr in atstrs[Z]:
                  outfile.write(ostr)

        outfile.close()

class veloc:
    """
    Class for handling velocities (stored in a.u.).
    """
    def read_file(self, file_path, file_type):
        """
        Read in the structure from a file.
        """
        self.file_path = file_path
        self.file_type = file_type

        infile = open(self.file_path,'r')
        line = infile.readline()
        vtmp = []

        if self.file_type == 'vtxyz': # read from tinker.dyn file
            inveloc = False
            while(line!=''):
                if 'Current Atomic Velocities' in line:
                    inveloc = True
                elif 'Accelerations' in line:
                    break
                elif inveloc:
                    words = [float(word.replace("D","E")) for word in line.split()]
                    vtmp.append([word/units.length['A']*units.time['fs']/1000 for word in words]) # convert to a.u. from Ang/ps
                line = infile.readline()

        elif self.file_type == 'vnx': # NX veloc file
           while(line!=''):
              vtmp.append([float(word) for word in line.split()])

              line = infile.readline()
        else:
           print 'type %s not supported for input'%self.file_type
           exit(1)

        infile.close()

        self.veloc = numpy.array(vtmp)
        #print self.veloc

    def read_struc(self, struc, scale=1.0):
        """
        Initialize with the coordinates of a structure file.
        scale - Optional scaling factor
        """
        self.veloc = struc.ret_3xN_matrix() / units.length['A'] * scale

    def write_veloc(self, file_path, file_type):
        wfile = open(file_path, 'w')

        if file_type == 'vnx':
            for at in self.veloc:
                wfile.write(" % 14.9f % 14.9f % 14.9f\n"%(at[0],at[1],at[2]))
        else:
           print 'type %s not supported for output'%file_type
           exit(1)

        wfile.close()