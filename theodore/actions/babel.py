from __future__ import print_function, division
import sys
from .. import theo_header, lib_struc
from .actions import Action


class Babel(Action):

    name = 'babel'

    _colt_description = 'Openbabel wrapper - conversion of coordinate files'

    _questions = """
    # Input molecular structure file
    infile = :: existing_file
    # Output molecular structure file
    outfile = :: file
    # File type for input file
    intype = :: str, optional, alias=i
    # File type for output file
    outtype = :: str, optional, alias=o
    """

    def run(infile, intype, outtype, outfile):
        theo_header.print_header(_colt_description)

        if intype in lib_struc.veloc_types: # special treatment of velocities
            veloc = struc_linalg.veloc()
            veloc.read_file(file_path=infile, file_type=intype)
            veloc.write_veloc(file_path=outfile,file_type=outtype)
        else:
            struc = lib_struc.structure()
            struc.read_file(file_path=infile, file_type=intype)
            struc.make_coord_file(file_path=outfile,file_type=outtype)
        
        print("Finished: file %s written."%outfile)
