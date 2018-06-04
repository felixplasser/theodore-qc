#!/usr/bin/env python2
"""
Plot arrows for dipole and quadrupole moments.
"""

import theo_header, units, lib_file, input_options, error_handler

class mom_options(input_options.write_options):
    def input(self):
        self.read_yn('Plot dipole moments', 'do_dip', True)
        if self['do_dip']:
            self.read_float('Scale factor for dipole moments', 'dip_scale', 2.0)
            self.read_float('Radius for dipole moments', 'dip_rad', 0.4)
        self.read_yn('Plot quadrupole moments', 'do_quad', True)
        if self['do_quad']:
            self.read_float('Scale factor for dipole moments', 'quad_scale', 1.0)
            self.read_float('Radius for quadrupole moments', 'quad_rad', 0.2)

    def write_afile(self, filen='arrows.vmd'):
        """
        File for write arrows for the different states.
        """
        af = open(filen, 'w')
        af.write(\
"""axes location Off
display projection Orthographic
display depthcue off
color Display Background white
menu graphics on
mol modstyle 0 0 Licorice 0.100000 30.000000 30.000000
""")

        sfile = lib_file.summ_file(self['ana_file'])
        ddict = sfile.ret_ddict()
        dfac = self['dip_scale']  * units.length['A']
        qfac = self['quad_scale'] * units.length['A']
        for state in sfile.ret_state_labels():
            sdict = ddict[state]

            # Dipole and quadrupole moments
            af.write('draw delete all\n')
            if self['do_dip']:
                if not 'mux' in sdict:
                    print " *** No dipole info found for state %s"%state
                else:
                    af.write('draw color green\n')
                    af.write('draw cylinder ')
                    self.vmd_coors(-.5 * dfac, .4 * dfac, sdict['mux'], sdict['muy'], sdict['muz'], af)
                    af.write('radius % .3f\n'%self['dip_rad'])

                    af.write('draw cone ')
                    self.vmd_coors( .4 * dfac, .6 * dfac, sdict['mux'], sdict['muy'], sdict['muz'], af)
                    af.write('radius % .3f\n'%(2*self['dip_rad']))

            if self['do_quad']:
                if not 'Qxx' in sdict:
                    print " *** No quadrupole info found for state %s"%state
                else:
                    tQxx = 2 * sdict['Qxx'] - sdict['Qyy'] - sdict['Qzz']
                    tQyy = 2 * sdict['Qyy'] - sdict['Qxx'] - sdict['Qzz']
                    tQzz = 2 * sdict['Qzz'] - sdict['Qyy'] - sdict['Qxx']
                    
                    if self.vmd_color(tQxx, af):
                        pQ = abs(tQxx)**.5 * units.length['A']
                        af.write('draw cylinder ')
                        self.vmd_coors(-.4 * dfac, .4 * dfac, pQ, 0., 0., af)
                        af.write('radius % .3f\n'%self['quad_rad'])
                        af.write('draw cone ')
                        self.vmd_coors(.4 * dfac, .6 * dfac, pQ, 0., 0., af)
                        af.write('radius % .3f\n'%(2*self['quad_rad']))
                        af.write('draw cone ')
                        self.vmd_coors(-.4 * dfac, -.6 * dfac, pQ, 0., 0., af)
                        af.write('radius % .3f\n'%(2*self['quad_rad']))
                    if self.vmd_color(tQyy, af):
                        pQ = abs(tQyy)**.5 * units.length['A']
                        af.write('draw cylinder ')
                        self.vmd_coors(-.4 * dfac, .4 * dfac, 0., pQ, 0., af)
                        af.write('radius % .3f\n'%self['quad_rad'])
                        af.write('draw cone ')
                        self.vmd_coors(.4 * dfac, .6 * dfac, 0., pQ, 0., af)
                        af.write('radius % .3f\n'%(2*self['quad_rad']))
                        af.write('draw cone ')
                        self.vmd_coors(-.4 * dfac, -.6 * dfac, 0., pQ, 0., af)
                        af.write('radius % .3f\n'%(2*self['quad_rad']))
                    if self.vmd_color(tQzz, af):
                        pQ = abs(tQzz)**.5 * units.length['A']
                        af.write('draw cylinder ')
                        self.vmd_coors(-.4 * dfac, .4 * dfac, 0., 0., pQ, af)
                        af.write('radius % .3f\n'%self['quad_rad'])
                        af.write('draw cone ')
                        self.vmd_coors(.4 * dfac, .6 * dfac, 0., 0., pQ, af)
                        af.write('radius % .3f\n'%(2*self['quad_rad']))
                        af.write('draw cone ')
                        self.vmd_coors(-.4 * dfac, -.6 * dfac, 0., 0., pQ, af)
                        af.write('radius % .3f\n'%(2*self['quad_rad']))                        

            af.write('render TachyonInternal state_%s.tga\n\n'%state)

        af.close()
        print "File %s written."%af.name

    def vmd_coors(self, fac1, fac2, x, y, z, af):
        af.write('{% .3f % .3f % .3f} '%(fac1 * x, fac1 * y, fac1 * z))
        af.write('{% .3f % .3f % .3f} '%(fac2 * x, fac2 * y, fac2 * z))

    def vmd_color(self, val, af, eps=1.e-3):
        if abs(val) < eps:
            return False
        elif val > 0.0:
            af.write('draw color blue\n')
            return True
        else:
            af.write('draw color red\n')
            return True            

if __name__=='__main__':
    import sys

    theo_header.print_header('Plotting of dipole and quadrupole moments')
    print 'draw_moments.py <tden_summ>'
    if len(sys.argv) < 2:
        raise error_handler.MsgError('Enter one argument')

    opt = mom_options('mom.in')
    opt['ana_file'] = sys.argv[1]
    opt.input()
    opt.write_afile()    

