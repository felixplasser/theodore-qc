#!/usr/bin/env python3
"""
Automatic plotting of vibrations with jmol.
"""
# Code adapted from jmol_MOs.py

from __future__ import print_function, division

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    lib_file = importer.lazy_import_as('..lib_file', 'lib_file')
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')


class jmol_vib_opts(input_options.write_options):
    def input(self):
        self.read_int('Index of first vibration to be plotted', 'st_ind', 7)
        self.read_int('Index of last vibration to be plotted',  'en_ind', 20)
        self.read_int('Width of vibration vectors', 'vwidth', 5)
        self.read_int('Scale of vibration vectors', 'vscale', 5)

        self.read_yn('Use "rotate best" command (only available in Jmol 14)', 'rot_best', True)
        self.read_yn('Additional custom rotation of the molecule?', 'rot_custom')
        if self['rot_custom']:
            self.read_float('Rotation around the x-axis', 'rot_x', 0.)
            self.read_float('Rotation around the y-axis', 'rot_y', 0.)
            self.read_float('Rotation around the z-axis', 'rot_z', 0.)
        self.read_str('Format of the output files (png, pngt, ...)', 'oformat', 'png')
        self.read_int('Width of images in output html file', 'width', 400)
        self.read_yn('Run Jmol?', 'run_jmol', False)

class vib_output:
    """
    Abstract base class for vib output.
    """
    def __init__(self, vibc, jopt):
        self.vibc = vibc
        self.jopt = jopt
        self.outstr = ''

    def output(self, ofileh):
        self.pre()
        self.print_vibs()
        self.post(ofileh)

    def vibpath(self, ivib, of='png'):
        if of == 'pngt': of = 'png'
        return "vib_%i.%s"%(ivib,of)

    def pre(self):
        raise error_handler.PureVirtualError()

    def print_vibs(self):
        raise error_handler.PureVirtualError()

    def post(self, ofileh):
        ofileh.write(self.outstr)

class vib_output_jmol(vib_output):
    """
    Vib output in standard jmol format.
    """
    def pre(self):
        self.outstr += '\nload "' + self.vibc.mldfile + '" FILTER "nosort"\n'
        if self.jopt['rot_best']:
            self.outstr += "rotate best\n"
        if self.jopt['rot_custom']:
            self.outstr += "rotate x %.3f\n"%self.jopt['rot_x']
            self.outstr += "rotate y %.3f\n"%self.jopt['rot_y']
            self.outstr += "rotate z %.3f\n"%self.jopt['rot_z']
        self.outstr += "background white\n" + "vector on\n"
        self.outstr += "vector %i\n"%self.jopt['vwidth']
        self.outstr += "vector scale %.2f\n"%self.jopt['vscale']

    def print_vibs(self):
        for ivib in self.vibc.viblist:
            self.outstr += "model %i\n"%(ivib+1)
            self.outstr += "write image %s \"%s\"\n"%(self.jopt['oformat'], self.vibc.vibpath(ivib, self.jopt['oformat']))

class vib_output_html(vib_output):
    """
    HTML file for visualizing vibrations created with jmol.
    """
    def pre(self):
        self.htable = lib_file.htmltable(ncol=3)

    def print_vibs(self):
        for ivib in self.vibc.viblist:
            el = '<img src="%s" "border="1" width="%i">'%(self.vibc.vibpath(ivib, self.jopt['oformat']), self.jopt['width'])
            el += '<br> Mode %i'%(ivib)
            self.htable.add_el(el)

    def post(self, ofileh):
        ofileh.write("<h2>Vibrations - %s</h2>\n"%self.vibc.mldfile)
        self.htable.close_table()
        ofileh.write(self.htable.ret_table())        

class vibcoll:
    """
    Output of molecular vibrations.
    """
    def __init__(self, st_ind, en_ind, mldfile):
        self.mldfile = mldfile
        self.viblist = list(range(st_ind, en_ind+1))

    def vibname(self, ivib):
        return str(ivib)

    def vibpath(self, ivib, of='png'):
        if of == 'pngt': of = 'png'
        return "vib_%i.%s"%(ivib,of)    


class JMolVibs(Action):

    name = 'jmol_vibs'
    
    _colt_description = 'Plotting of vibrations in Jmol'

    _user_input = """
    # Molden file with info about vibrations
    vibfile = :: existing_file
    """

    _lazy_imports = LazyImporter({
            '..error_handler': 'error_handler',
            '..lib_file': 'lib_file',
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
    })

    def run(vibfile):

        theo_header.print_header(__class__._colt_description)

        jopt = jmol_vib_opts('jmol.in')
        jopt.input()
        
        jo = lib_file.wfile('jmol_vibs.spt')
        ho = lib_file.htmlfile('vibs.html')
        
        ho.pre('Vibrations')
        
        vibc = vibcoll(jopt['st_ind'], jopt['en_ind'], vibfile)

        vibout = vib_output_jmol(vibc, jopt)
        vibout.output(jo)
        
        vibh = vib_output_html(vibc, jopt)
        vibh.output(ho)

        ho.post(lvprt=1)
        jo.post(lvprt=1)
        if jopt['run_jmol']:
            import subprocess
            print("Running jmol ...")

            subprocess.call(["jmol", "-n", jo.name])
        else:
            print("  -> Now simply run \"jmol -n %s\" to plot all the orbitals.\n"%jo)
