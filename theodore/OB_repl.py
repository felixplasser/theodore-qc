"""
version 1.0.1
author: Felix Plasser, University of Vienna, Institute for Theoretical Chemistry
Waehringerstr. 17, 1090, Vienna, Austria
usage: OB_repl is intended as a replacement that supplies some of the functionalities of the python-openbabel
   package in case this is not available.
"""

from __future__ import print_function, division

import sys
from . import units

Z_mass_dict = {1:1.008,6:12.011,7:14.007,8:16.,15:30.97,16:32.066,17:35.453,77:192.22}
Z_exact_mass_dict = {1:1.00782504,6:12.,7:14.00307401,8:15.99491464,17:34.96885273,77:192.962917}
from .atominfo import symbol_Z_dict, Z_symbol_dict

class OBConversion:
    """
    Emulate the openbabel OBConversion class with support for a reduced number of file formats.
    """
    def __init__(self):
        self.in_format = None
        self.out_format = None

    def SetInFormat(self,iformat):
        self.in_format = iformat

        return True

    def SetOutFormat(self,oformat):
        self.out_format = oformat

        return True

    def ReadFile(self,mol,file):
        """
        Read information from a file into mol.
        """
        if self.in_format == 'xyz':
            self.read_xyz(mol, file)
        elif self.in_format == 'tmol':
            self.read_tmol(mol, file)
        elif self.in_format == None:
            print("Input format for %s not set!"%file)
            sys.exit(14)
        else:
            print("File format %s not supported for input!"%self.in_format)
            print("Install python-openbabel for complete support.")
            sys.exit(15)

        return True

    def read_xyz(self,mol,file):
        """
        Read file in xyz format (Angstrom).
        """
        infile = open(file,'r')

        # the first two lines do not contain any coordinate information
        line = infile.readline()
        line = infile.readline()
        line = infile.readline()


        while(line!=''):
              words = line.split()
              obatom = OBAtom()
              obatom.SetAtomicNum(symbol_Z_dict[words[0]])
              coords = [float(word) for word in words[1:4]]
              obatom.SetVector(*coords)

              mol.AddAtom(obatom)
              line = infile.readline()

    def read_tmol(self, mol, file):
        """
        Read Turbomole format (Bohr)
        """
        infile = open(file,'r')

        line = infile.readline()
        line = infile.readline()

        while(not '$' in line):
              words = line.split()
              obatom = OBAtom()
              obatom.SetAtomicNum(symbol_Z_dict[words[3][0].upper()+words[3][1:]])
              coords = [float(word) * units.length['A'] for word in words[0:3]]
              obatom.SetVector(*coords)

              mol.AddAtom(obatom)
              line = infile.readline()

    def WriteFile(self,mol,file):
        """
        Print structure to a file.
        """
        if self.out_format == 'xyz':
            self.write_xyz(mol,file)
        elif self.out_format == 'tmol':
            self.write_tmol(mol,file)
        elif self.out_format == None:
            print("Output format for %s not set!"%file)
            sys.exit(14)
        else:
            print("File format %s not supported for output!"%self.out_format)
            print("Install python-openbabel for complete support.")
            sys.exit(15)

    def write_xyz(self,mol,file):
        outfile = open(file,'w')
        num_at = mol.NumAtoms()

        outfile.write('%i\n\n'%num_at)

        for ind in range(1, num_at+1):
            obatom  = mol.GetAtom(ind)
            outstr  = '%2s'%Z_symbol_dict[obatom.GetAtomicNum()]
            outstr += '% 14.8f'%(obatom.x())
            outstr += '% 14.8f'%(obatom.y())
            outstr += '% 14.8f'%(obatom.z())
            outfile.write(outstr+'\n')

        outfile.close()

    def write_tmol(self,mol,file):
        outfile = open(file,'w')
        num_at = mol.NumAtoms()

        outfile.write('$coord\n')

        for ind in range(1, num_at+1):
            obatom  = mol.GetAtom(ind)

            outstr  = '% 14.8f'%(obatom.x() / units.length['A'])
            outstr += '% 14.8f'%(obatom.y() / units.length['A'])
            outstr += '% 14.8f'%(obatom.z() / units.length['A'])
            outstr += '%2s'%Z_symbol_dict[obatom.GetAtomicNum()]
            outfile.write(outstr+'\n')

        outfile.write('$end\n')
        outfile.close()

class OBMol:
    def __init__(self):
        self.atoms = []

    def AddAtom(self,Atom):
        self.atoms.append(Atom)

    def GetAtom(self,ind):
        """
        Return an atom, indices start with 1.
        """
        return self.atoms[ind-1]

    def NumAtoms(self):
        return len(self.atoms)

class OBAtom:
    def __init__(self):
        self.AtomicNum = None
        self.Vector = [None,None,None]

    def SetAtomicNum(self,AtomicNum):
        self.AtomicNum = AtomicNum

    def GetAtomicNum(self):
        return self.AtomicNum

    def SetVector(self,*Vector):
        self.Vector=Vector

    def GetVector(self):
        return self.Vector

    def x(self):
        return self.Vector[0]

    def y(self):
        return self.Vector[1]

    def z(self):
        return self.Vector[2]

    def GetExactMass(self):
        return Z_exact_mass_dict[self.AtomicNum]

class OBAtomAtomIter:
    def __init__(self, atom):
        print("\n Functionality not supported!")
        print("Install python-openbabel for complete support.")
        sys.exit(16)
