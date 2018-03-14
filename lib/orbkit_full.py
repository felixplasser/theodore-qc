"""
This is an interface to the orbkit, an external post-processing toolbox.
https://orbkit.github.io/
orbkit should be part of this distribution. Otherwise,
download and install orbkit if you want to use the functions.
Author: Gunter Hermann
"""

import dens_ana_base
import numpy,tempfile

# Import orbkit modules
from orbkit import read,grid,core,output,options
from orbkit.qcinfo import QCinfo
from orbkit.core import l_deg,lquant
from orbkit.display import display
from orbkit.detci import ci_core

# Disable orbkit terminal output for each run
options.quiet = True
options.no_log = True

class lib_orbkit:

    def __init__(self):
        self.slice_length = 1e4

    def orbkit_geo_ao_conversion(self,mos):

        # Read geader with molden writer from orbkit
        tmp = tempfile.NamedTemporaryFile()

        # Open the file for writing.
        with open(tmp.name, 'w') as f:
            f.write(mos.header)
            f.write('[mo]\nSym= 1a\nEne= 0.0\nSpin= Alpha\nOccup= 2.0\n')
        qc = read.read_molden(tmp.name)

        return qc

    def orbkit_nto_conversion(self,U,lam,Vt,mos,qc,minlam=1e-3):

        # MO conversion
        qc.mo_spec = []
        Csub = mos.mo_mat.transpose()[:len(U)]
        U_mat_t =  numpy.dot(Csub.transpose(), U)
        U_mat_t = (U_mat_t).transpose()
        V_mat_t = numpy.dot(Vt,mos.mo_mat.transpose())
        UV_t = [iU for iU in reversed(U_mat_t)] + [iV for iV in V_mat_t]
        lam2 = [-vlam for vlam in reversed(lam)] + [vlam for vlam in lam] + [0.] * (len(Vt) - len(lam))
        LUMO = (numpy.array(mos.occs) > 0).sum()
        for i in range(len(lam2)):
            if abs(lam2[i]) > minlam:
                qc.mo_spec.append({'coeffs': UV_t[i],
                                   'energy': lam2[i],
                                   'occ_num': lam2[i],
                                   'sym': 'NTO_h' if i < LUMO else 'NTO_p',
                                   'spin': 'Alpha'
                                 })
        return qc

    def orbkit_mo_conversion(self,mos,qc):

        # MO conversion
        qc.mo_spec = []
        mo_mat = mos.mo_mat.T
        for i in range(len(mo_mat)):
            try:
                symtmp = mos.syms[i]
            except IndexError:
                symtmp = 'X'
            qc.mo_spec.append({'coeffs': mo_mat[i],
                                'energy': mos.ens[i],
                                'occ_num': mos.occs[i],
                                'sym': symtmp,
                                'spin': 'Alpha'
                              })
        return qc

    def orbkit_grid(self,qc):

        # Initialize grid
        grid.adjust_to_geo(qc,extend=5.0,step=0.4)
        grid.grid_init(force=True)
        self.slice_length = grid.N_[1]*grid.N_[2]/2

    def compute_MOs(self,qc,numproc=4):

        self.orbkit_grid(qc)
        molist = core.rho_compute(qc,calc_mo=True,slice_length=self.slice_length,drv=None,numproc=numproc)

        return molist

    def compute_p_h_dens(self,state, U, lam, Vt, mos, minlam=1e-3,numproc=4):

        print("Calculating particle/hole density with orbkit for state %s" % (state['name']))

        # Data conversion from TheoDORE to orbkit
        qc = self.orbkit_geo_ao_conversion(mos)
        qc = self.orbkit_nto_conversion(U,lam,Vt,mos,qc,minlam=minlam)

        # Calculate MOs
        molist = self.compute_MOs(qc,numproc=numproc)

        # Calculate hole and particle density
        rho_p = numpy.zeros(tuple(grid.N_))
        rho_h = numpy.zeros(tuple(grid.N_))
        for i in range(len(qc.mo_spec)):
            if 'p' in qc.mo_spec[i]['sym']:
                rho_p += abs(qc.mo_spec[i]['energy'])*molist[i]*molist[i]
            else:
                rho_h += abs(qc.mo_spec[i]['energy'])*molist[i]*molist[i]

        # Reshape particle and hole density and write cube-files
        cube_ids = []
        fid = 'rho_p_%s'% state['name'].replace('(', '-').replace(')', '-')
        cube_ids.append(fid)
        output.cube_creator(rho_p,fid,qc.geo_info,qc.geo_spec)
        fid = 'rho_h_%s' % state['name'].replace('(', '-').replace(')', '-')
        output.cube_creator(rho_h,fid,qc.geo_info,qc.geo_spec)
        cube_ids.append(fid)

        return cube_ids

    def compute_rho_0_n(self,state_list,mos,numproc=4):

        # Data conversion from TheoDORE to orbkit
        qc = self.orbkit_geo_ao_conversion(mos)
        qc = self.orbkit_mo_conversion(mos,qc)

        # Calculate MOs
        molist = self.compute_MOs(qc)

        zero = [[],[]]
        cube_fids = []
        for i in range(len(state_list)):
            sing = [[],[]]
            print("Transition density between ground state and excited state %s" % (state_list[i]['name']))
            for j in range(state_list[i]['tden'].shape[0]):
                for k in range(state_list[i]['tden'].shape[1]):
                  if abs(state_list[i]['tden'][j,k]) >= 1e-8:
                    sing[0].append(state_list[i]['tden'][j,k])
                    sing[1].append([j,k])
            rho0n = ci_core.rho(zero,sing,molist,slice_length=self.slice_length,numproc=numproc)
            fid = 'rho_0_%s' % (state_list[i]['name'].replace('(', '-').replace(')', '-'))
            output.cube_creator(rho0n,fid,qc.geo_info,qc.geo_spec)
            cube_fids.append(fid)

        return cube_fids

    def cube_file_creator(self,state, U, lam, Vt, mos, minlam=1e-3, numproc=4):

        print("Calculating NTOs as cube files with orbkit for state %s" % (state['name']))

        # Data conversion from TheoDORE to orbkit
        qc = self.orbkit_geo_ao_conversion(mos)
        qc = self.orbkit_nto_conversion(U,lam,Vt,mos,qc,minlam=minlam)
        # Calculate MOs
        molist = self.compute_MOs(qc)

        cube_ids = []
        # Create Cube-files
        for i in range(len(molist)):
            fid = '%s_%s_%.2f'% (qc.mo_spec[i]['sym'], state['name'].replace('(', '-').replace(')', '-'),abs(qc.mo_spec[i]['occ_num']))
            cube_ids.append(fid)
            output.cube_creator(molist[i],fid,qc.geo_info,qc.geo_spec)
        return cube_ids

    def vmd_network_creator(self,filename='',cube_ids=[],isovalue=0.01):
        cube_ids = ['%s.cb' % i for i in cube_ids]
        output.vmd_network_creator(filename,cube_files=cube_ids,render=False,iso=(-isovalue,isovalue))
