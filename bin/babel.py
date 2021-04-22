from __future__ import print_function, division
import sys
from .. import theo_header, lib_struc


def babel():
    theo_header.print_header('Openbabel wrapper - conversion of coordinate files')
    
    print("""\
    Usage: molecular structure conversion based on the openbabel package
         - increased support for Columbus (col) and Tinker (txyz2) formats
         - velocity conversion Newton-X (vnx) and Tinker (vtxyz)
    
    Syntax: babel.py <intype> <infile> <outtype> <outfile>
       e.g. babel.py tmol coord xyz coord.xyz
       or   babel.py <infile> <outfile>
    
    """)
    
    if len(sys.argv) == 4+1:
        (intype,infile,outtype,outfile) = sys.argv[1:]
    elif len(sys.argv) == 2+1:
        (infile,outfile) = sys.argv[1:]
        intype = outtype = None
    else:
        print("""  Supported file types:
        col -- Columbus and Newton-X format
        colr -- Columbus format, atoms reordered
        txyz2 -- Tinker format with reading possibility (verify the atom types in the output)
        vnx -- veloctiy (Newton-X format)
        ntxyz -- velocity (Tinker format)
      additionally all formats from openbabel are included
        type 'babel -H' for a complete list""")
    
        print('\nTwo or four arguments required.')
        sys.exit()
    
    if intype in lib_struc.veloc_types: # special treatment of velocities
        veloc = struc_linalg.veloc()
        veloc.read_file(file_path=infile, file_type=intype)
        veloc.write_veloc(file_path=outfile,file_type=outtype)
    else:
        struc = lib_struc.structure()
        struc.read_file(file_path=infile, file_type=intype)
        struc.make_coord_file(file_path=outfile,file_type=outtype)
    
    print("Finished: file %s written."%outfile)
