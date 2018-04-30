#!/usr/bin/env python2
"""
Automatic plotting of densities or MOs with vmd.
"""

import error_handler, theo_header, input_options, lib_struc, lib_file

class vmd_options(input_options.write_options):
    def vmd_input(self):
        self.read_int('How many isovalues (1 or 2)?', 'niso', 2)
        self.read_float('First isovalue', 'iso1', 0.003)
        if self['niso'] >= 2:
            self.read_float('Second isovalue', 'iso2', self['iso1']/3.)

        if self.ret_yn('Adjust detailed options?', False):
            self.read_str('Name of the file used to load data into VMD', 'lfile', 'load_all.vmd')
            self.read_str('Name of the file used to plot in VMD', 'pfile', 'plot_all.vmd')
            self.read_str('Name of the file used to call GIMP convert', 'cfile', 'convert.bash')
            self.read_str('Name of the HTML file with the plots', 'hfile', 'vmd_plots.html')
            self.read_int('Width of images in output html file', 'width', 400)
        else:
            self['lfile'] = 'load_all.vmd'
            self['pfile'] = 'plot_all.vmd'
            self['cfile'] = 'convert.bash'
            self['hfile'] = 'vmd_plots.html'
            self['width'] = 400

    def write_lfile(self, pltfiles):
        """
        File for loading data.
        """
        lf = open(self['lfile'], 'w')
        lf.write("""material change opacity Glass3 0.150000
material change diffuse Glass3 0.10000
axes location Off
display projection Orthographic
display rendermode GLSL
display depthcue off
color Display Background white
menu graphics on
material change diffuse Ghost 0.000000
material change ambient Ghost 0.300000
material change opacity Ghost 0.100000
material change shininess Ghost 0.000000
mol modstyle 0 0 Licorice 0.100000 30.000000 30.000000
""")
        lf.write("""mol addrep 0
mol addrep 0
mol modmaterial 1 0 AOShiny
mol modmaterial 2 0 AOShiny
mol modstyle 1 0 Isosurface  %.5f 0 0 0 1 1
mol modstyle 2 0 Isosurface -%.5f 0 0 0 1 1
mol modcolor 1 0 ColorID 0
mol modcolor 2 0 ColorID 1
"""%(self['iso1'], self['iso1']))

        if self['niso'] >= 2:
            lf.write("""mol addrep 0
mol addrep 0
mol modmaterial 3 0 Glass3
mol modmaterial 4 0 Glass3
mol modstyle 3 0 Isosurface  %.5f 0 0 0 1 1
mol modstyle 4 0 Isosurface -%.5f 0 0 0 1 1
mol modcolor 3 0 ColorID 0
mol modcolor 4 0 ColorID 1
"""%(self['iso2'], self['iso2']))

        struc = lib_struc.structure()
        for pltf in pltfiles:
            ftyp = struc.guess_file_type(pltf)
            lf.write("mol addfile %s type %s\n"%(pltf,ftyp))

        lf.close()
        print "File %s written."%lf.name

    def write_pfile(self, pltfiles):
        """
        File used for plotting.
        """
        pf = open(self['pfile'], 'w')
        for iplt, pltf in enumerate(pltfiles):
            pf.write("mol modstyle 1 0 Isosurface  %.5f %i 0 0 1 1\n"%(self['iso1'], iplt))
            pf.write("mol modstyle 2 0 Isosurface -%.5f %i 0 0 1 1\n"%(self['iso1'], iplt))
            if self['niso'] >= 2:
                pf.write("mol modstyle 3 0 Isosurface  %.5f %i 0 0 1 1\n"%(self['iso2'], iplt))
                pf.write("mol modstyle 4 0 Isosurface -%.5f %i 0 0 1 1\n"%(self['iso2'], iplt))
            pf.write("render TachyonInternal %s.tga\n"%pltf)

        pf.close()
        print "File %s written."%pf.name

    def write_cfile(self, pltfiles):
        """
        File for file conversion.
        """
        cf = open(self['cfile'], 'w')

        cf.write('#!/bin/bash\n')
        for pltf in pltfiles:
            cf.write("convert %s.tga %s.png\n"%(pltf, pltf))
            cf.write("rm %s.tga\n"%pltf)

        cf.close()
        print "File %s written."%cf.name

    def write_hfile(self, pltfiles):
        """
        HTML File.
        """
        ho = lib_file.htmlfile(self['hfile'])
        ho.pre('VMD plots')

        ht = lib_file.htmltable(ncol=4)
        for pltf in pltfiles:
            el = '<img src="%s.png" "border="1" width="%i">'%(pltf, self['width'])
            el += '<br> %s'%pltf
            ht.add_el(el)

        ht.close_table()
        ho.write(ht.ret_table())

        ho.post(lvprt=1)

def run():
    print 'vmd_plots.py [<pltfile> [<pltfile2> ...]]\n'

    pltfiles = sys.argv[1:]

    if len(pltfiles) == 0:
        raise error_handler.MsgError('No file specified')

    vopt = vmd_options('vmd.in')
    vopt.vmd_input()
    vopt.write_lfile(pltfiles)
    vopt.write_pfile(pltfiles)
    vopt.write_cfile(pltfiles)
    vopt.write_hfile(pltfiles)

    print "Converting coordinate file ..."
    struc = lib_struc.structure()
    try:
        struc.read_file(file_path=pltfiles[0], file_type=None)
        struc.make_coord_file(file_path='coord.xyz',file_type='xyz',lvprt=1)
    except:
        print "*** WARNING: The coordinate file coord.xyz could not be created. ***"
        print   "    Please create this file yourself.\n\n"

    print """
Files created. Now do the following:
1. vmd coord.xyz
2.   File - Load Visualization State - %s
3.   Adjust the perspective
4.   File - Load Visualization State - %s
5. bash %s
6. Open in browser: %s
"""%(vopt['lfile'], vopt['pfile'], vopt['cfile'], vopt['hfile'])

if __name__=='__main__':
    import sys

    theo_header.print_header('Automatic plotting in VMD')
    run()
