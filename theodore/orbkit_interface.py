"""
This is an interface to the orbkit, an external post-processing toolbox.
https://orbkit.github.io/
Download and install orbkit if you want to use the functions.
"""
from __future__ import print_function, division

orbkit_avail = True
try:
    import orbkit
except ImportError:
    orbkit_avail = False
    print(" *** Orbkit installation not found. ***")
    print(" Install orbkit for extended plotting capabilities.\n")

if orbkit_avail:
    from .orbkit_full import lib_orbkit
else:
    class lib_orbkit:
        """
        This is a fake class just to make sure that all the calls are defined.
        """
        def __init__(self):
            print("\n*** WARNING: orbkit not available! ***\n    Plotting not possible as specified.\n")
        def orbkit_geo_ao_conversion(self,mos):
            pass
        def orbkit_nto_conversion(self,U,lam,Vt,mos,qc,minlam=1e-3):
            pass
        def orbkit_mo_conversion(self,mos,qc):
            pass
        def orbkit_grid(self,qc):
            pass
        def compute_MOs(self,qc,numproc=4):
            pass
        def compute_p_h_dens(self,state, U, lam, Vt, mos, minlam=1e-3,numproc=4):
            pass
        def compute_rho_0_n(self,state_list,mos,numproc=4):
            pass
        def cube_file_creator(self,state, U, lam, Vt, mos, minlam=1e-3, numproc=4):
            pass
        def vmd_network_creator(self,filename='',cube_ids=[],isovalue=0.01):
            pass
