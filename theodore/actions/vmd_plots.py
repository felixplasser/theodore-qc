"""
Automatic plotting of densities or MOs with vmd.
"""
from __future__ import print_function, division
import os
import sys
from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    lib_struc = importer.lazy_import_as('..lib_struc', 'lib_struc')
    lib_file = importer.lazy_import_as('..lib_file', 'lib_file')
    lib_util = importer.lazy_import_as('..lib_util', 'lib_util')


class vmd_options(input_options.write_options):
    def vmd_input(self):
        self.read_yn('Compute volume integrals over cube files for isovalues?', 'do_vol', False)
        self.read_yn('Use special DNTO mode?', 'dnto', False)
        if self['dnto']:
            self['niso'] = 2
            if not self['do_vol']:
                self.read_float('Isovalue for conditional density', 'iso1', 0.003)
                self.read_float('Isovalue for probe density', 'iso2', self['iso1'])
            else:
                self.read_float('Volume integral for conditional density', 'iso1', 0.75)
                self.read_float('Volume integral for probe density', 'iso2', self['iso1'])
            self.read_str('VMD Material for conditional density [AOShiny, EdgyGlass, ...]', 'mat1', 'AOShiny')
            self.read_str('VMD Material for probe density [Glass1, AOEdgy, ...]', 'mat2', 'Glass1')
        else:
            self.read_int('How many isovalues (1 or 2)?', 'niso', 2)
            if not self['do_vol']:
                self.read_float('First isovalue', 'iso1', 0.003)
            else:
                self.read_float('Volume integral for first density', 'iso1', 0.5)
            self.read_str('VMD Material for first density', 'mat1', 'AOShiny')
            if self['niso'] >= 2:
                if not self['do_vol']:
                    self.read_float('Second isovalue', 'iso2', self['iso1']/3.)
                else:
                    self.read_float('Volume integral for second density', 'iso2', 0.9)
                self.read_str('VMD Material for second density', 'mat2', 'Glass3')
            else:
                self['iso2'] = 100.

        self.read_int('Width of images in output html file', 'width', 400)
        self.read_int('Number of columns in the output html file', 'ncol', 4)
        if self.ret_yn('Adjust file names?', False):
            self.read_str('Name of the file used to load data into VMD', 'lfile', 'load_all.vmd')
            self.read_str('Name of the file used to plot in VMD', 'pfile', 'plot_all.vmd')
            self.read_str('Name of the file used to call GIMP convert', 'cfile', 'convert.bash')
            self.read_str('Name of the HTML file with the plots', 'hfile', 'vmd_plots.html')
        else:
            self['lfile'] = 'load_all.vmd'
            self['pfile'] = 'plot_all.vmd'
            self['cfile'] = 'convert.bash'
            self['hfile'] = 'vmd_plots.html'

    def mod_pltfiles(self, pltfiles):
        """
        Separate plotfiles for DNTO mode.
        """
        hfiles = []
        auxhfiles = []
        efiles = []
        auxefiles = []
        for pltf in pltfiles:
            if 'rho_p' in pltf and not 'elec' in pltf:
                hfiles.append(pltf)
                auxhfiles.append(pltf.replace('rho_p', 'rho_h'))
            elif 'rho_h' in pltf and not 'hole' in pltf:
                efiles.append(pltf)
                auxefiles.append(pltf.replace('rho_h', 'rho_p'))
        return hfiles + efiles, auxhfiles + auxefiles

    def write_lfile(self, pltfiles, auxfiles=[]):
        """
        File for loading data.
        """
        lf = open(self['lfile'], 'w')
        lf.write("""material change opacity Glass3 0.150000
material change diffuse Glass3 0.10000
axes location Off
display projection Orthographic
display depthcue off
color Display Background white
menu graphics on
mol modstyle 0 0 Licorice 0.100000 30.000000 30.000000
""")
# material change diffuse Ghost 0.000000
# material change ambient Ghost 0.300000
# material change opacity Ghost 0.100000
# material change shininess Ghost 0.000000

        iso1 = 0.001 if self['do_vol'] else self['iso1']
        lf.write("""mol addrep 0
mol addrep 0
mol modmaterial 1 0 %s
mol modmaterial 2 0 %s
mol modstyle 1 0 Isosurface  %.5f 0 0 0 1 1
mol modstyle 2 0 Isosurface -%.5f 0 0 0 1 1
mol modcolor 1 0 ColorID 0
mol modcolor 2 0 ColorID 1
"""%(self['mat1'], self['mat1'], iso1, iso1))

        if self['niso'] >= 2:
            iso2 = 0.001 if self['do_vol'] else self['iso2']
            lf.write("""mol addrep 0
mol addrep 0
mol modmaterial 3 0 %s
mol modmaterial 4 0 %s
mol modstyle 3 0 Isosurface  %.5f 0 0 0 1 1
mol modstyle 4 0 Isosurface -%.5f 0 0 0 1 1
mol modcolor 3 0 ColorID 0
mol modcolor 4 0 ColorID 1
"""%(self['mat2'], self['mat2'], iso2, iso2))

        struc = lib_struc.structure()
        for pltf in pltfiles:
            ftyp = struc.guess_file_type(pltf)
            lf.write("mol addfile %s type %s\n"%(pltf,ftyp))
        for pltf in auxfiles:
            ftyp = struc.guess_file_type(pltf)
            lf.write("mol addfile %s type %s\n"%(pltf,ftyp))

        lf.close()
        print("File %s written."%lf.name)

    def write_pfile(self, pltfiles, auxfiles=[]):
        """
        File used for plotting.
        """
        iso1 = self['iso1']
        iso2 = self['iso2']
        pf = open(self['pfile'], 'w')
        for iplt, pltf in enumerate(pltfiles):
            if self['do_vol']:
                isovals = lib_util.cube_file(pltf).ret_isovals([self['iso1'], self['iso2']], lvprt=1)
                iso1 = isovals[0]
            pf.write("mol modstyle 1 0 Isosurface  %.5f %i 0 0 1 1\n"%(iso1, iplt))
            pf.write("mol modstyle 2 0 Isosurface -%.5f %i 0 0 1 1\n"%(iso1, iplt))
            if self['dnto']:
                if self['do_vol']:
                    iso2 = lib_util.cube_file(auxfiles[iplt]).ret_isovals([self['iso2']], lvprt=1)[0]
                pf.write("mol modstyle 3 0 Isosurface  %.5f %i 0 0 1 1\n"%(iso2, iplt + len(pltfiles)))
                pf.write("mol modstyle 4 0 Isosurface -%.5f %i 0 0 1 1\n"%(iso2, iplt + len(pltfiles)))
            elif self['niso'] >= 2:
                if self['do_vol']:
                    iso2 = isovals[1]
                pf.write("mol modstyle 3 0 Isosurface  %.5f %i 0 0 1 1\n"%(iso2, iplt))
                pf.write("mol modstyle 4 0 Isosurface -%.5f %i 0 0 1 1\n"%(iso2, iplt))
            pf.write("render TachyonInternal %s.tga\n"%pltf)

        pf.close()
        print("File %s written."%pf.name)

    def write_cfile(self, pltfiles):
        """
        File for file conversion.
        """
        cf = open(self['cfile'], 'w')

        cf.write('#!/bin/bash\n')
        for pltf in pltfiles:
            cf.write("convert %s.tga %s.png && "%(pltf, pltf))
            cf.write("rm %s.tga\n"%pltf)

        cf.close()
        print("File %s written."%cf.name)

    def write_hfile(self, pltfiles):
        """
        HTML File.
        """
        ho = lib_file.htmlfile(self['hfile'])
        ho.pre('VMD plots')

        ht = lib_file.htmltable(ncol=self['ncol'])
        for pltf in pltfiles:
            el = '<img src="%s.png" "border="1" width="%i">'%(pltf, self['width'])
            el += '<br> %s'%pltf
            ht.add_el(el)

        ht.close_table()
        ho.write(ht.ret_table())

        ho.post(lvprt=1)

class VMDPlots(Action):

    name = 'vmd_plots'

    _colt_description = 'Automatic plotting of cube files in VMD'

    _user_input = """
    # List of cube files (or other format VMD can read)
    pltfiles = :: list(existing_file)
    """

    _lazy_imports = LazyImporter({
            '..error_handler': 'error_handler',
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
            '..lib_struc': 'lib_struc',
            '..lib_file': 'lib_file',
            '..lib_util': 'lib_util',
    })

    def run(pltfiles):
        theo_header.print_header(title=__class__._colt_description)

        print('%i Files analyzed:' % len(pltfiles), end=' ')
        print(", ".join(os.path.basename(filename) for filename in pltfiles))

        vopt = vmd_options('vmd.in')
        vopt.vmd_input()
        auxfiles = []
        if vopt['dnto']:
            pltfiles, auxfiles = vopt.mod_pltfiles(pltfiles)

        vopt.write_lfile(pltfiles, auxfiles)
        vopt.write_pfile(pltfiles, auxfiles)
        vopt.write_cfile(pltfiles)
        vopt.write_hfile(pltfiles)

        print("Converting coordinate file ...")
        struc = lib_struc.structure()
        try:
            struc.read_file(file_path=pltfiles[0], file_type=None)
            struc.make_coord_file(file_path='coord.xyz',file_type='xyz',lvprt=1)
        except:
            print("*** WARNING: The coordinate file coord.xyz could not be created. ***")
            print("    Please create this file yourself.\n\n")

        print("""
Files created. Now do the following:
1. vmd coord.xyz
2.   File - Load Visualization State - %s
3.   Adjust the perspective
4.   File - Load Visualization State - %s
5. bash %s
6. Open in browser: %s
"""%(vopt['lfile'], vopt['pfile'], vopt['cfile'], vopt['hfile']))
