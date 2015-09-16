#!/usr/bin/python
"""
Automatic plotting of MOs with jmol.
"""
# TODO: input options

import lib_mo, error_handler, lib_file, theo_header, input_options

class mo_output:
    """
    Abstract base class for MO output.
    """
    def __init__(self, moc, jopt):
        self.moc = moc
        self.jopt = jopt
        self.outstr = ''
        
    def output(self, ofileh):
        self.pre()
        self.print_mos()
        self.post(ofileh)
        
    def mopath(self, mo):
        return "MO_%s.png"%mo
    
    def pre(self):
        raise error_handler.PureVirtualError()

    def print_mos(self):
        raise error_handler.PureVirtualError()

    def post(self, ofileh):
        ofileh.write(self.outstr)
        
class mo_output_jmol(mo_output):
    """
    MO output in standard jmol format.
    """
    def pre(self):
        if self.moc.mldfile != "":
            self.outstr += '\nload "' + self.moc.mldfile + '" FILTER "nosort"\n'
        self.outstr += "mo titleformat ''\n"
        if self.jopt['rot_best']:
            self.outstr += "rotate best\n"
        if self.jopt['rot_custom']:
            self.outstr += "rotate x %.3f\n"%self.jopt['rot_x']
            self.outstr += "rotate y %.3f\n"%self.jopt['rot_y']
            self.outstr += "rotate z %.3f\n"%self.jopt['rot_z']
        self.outstr += "background white\n" + "mo fill\n"
        self.outstr += "mo cutoff %.3f\n\n"%self.jopt['cutoff']
        
    def print_mos(self):
        for imo in self.moc.molist:
            self.outstr += "mo %s\n"%self.moc.moname(imo)
            self.outstr += "write image png \"%s\"\n"%(self.moc.mopath(imo))            
        
class mo_output_html(mo_output):
    """
    HTML file for visualizing the MOs created with jmol.
    """
    def pre(self):        
        self.htable = lib_file.htmltable(ncol=2)
        
    def print_mos(self):
        for imo in self.moc.molist:            
            el = '<img src="%s" "border="1" width="%i">'%(self.moc.mopath(imo), self.jopt['width'])
            el += self.moc.mo_extra(imo,pref="<br> MO %i:"%imo)
            self.htable.add_el(el)
        
    def post(self, ofileh):
        ofileh.write("<h2>Orbitals - %s</h2>\n"%self.moc.mldfile)
        self.htable.close_table()
        ofileh.write(self.htable.ret_table())        
        
class mo_output_tex(mo_output):
    """
    tex file for visualizing the MOs created with jmol.
    """
    def pre(self):
        self.ltable = lib_file.latextable(ncol=2)

    def print_mos(self):
        moex = []
        for imo in self.moc.molist:
            el = "\\incMO{%s}"%(self.moc.mopath(imo))
            moex += [self.moc.mo_extra(imo)]
            lastcol = self.ltable.add_el(el)
            if lastcol:
                self.ltable.add_row(moex)
                moex = []
            
    def post(self, ofileh):
        self.ltable.close_table()
        ofileh.write(self.ltable.ret_table()) 
        
class mocoll:
    def __init__(self, st_ind, en_ind, mldfile=""):
        self.mldfile = mldfile
        
        if not mldfile=="":
            self.moset = lib_mo.MO_set_molden(mldfile)
            self.moset.read(lvprt=0)
            maxmo = self.moset.ret_num_mo()
        else:
            maxmo = 10000
              
        self.molist = range(st_ind-1, min(en_ind, maxmo))
        
    def moname(self, imo):
        return str(imo+1)
        
    def mopath(self, imo):
        return "%sMO_%s.png"%(self.ret_label(),self.moname(imo))
        
    def mo_extra(self, imo, pref="", postf=""):
        if self.mldfile=="":
            return ""
        else:
            sym = self.moset.ret_sym(imo)
            try:
                ene, occ = self.moset.ret_eo(imo)
            except:
                print "\n ERROR: imo = %i"%imo
                print "ens:", self.moset.ens
                print "occs:", self.moset.occs
                raise
            return "%s %5s %.4f / %.4f %s"%(pref,sym,ene,occ,postf)
            
    def ret_label(self):
        if self.mldfile=="":
            return ""
        else:
            #return self.mldfile.split('.')[0] + '_'
            return self.mldfile.replace('.','-') + '_'
        
class mocollf(mocoll):
    """
    For frontier MOs.
    """
    def __init__(self, en_ind, mldfile=""):
        mocoll.__init__(self, 1, 2*en_ind, mldfile)
        
    def moname(self, imo):
        if imo==0:
            return "homo"
        elif imo==1:
            return "lumo"
        elif imo%2==0:
            return "homo-%i"%(imo/2)
        else:
            return "lumo+%i"%(imo/2)
            
    def mo_extra(self, imo, pref="", postf=""):
        if self.mldfile=="":
            return ""
        else:
            ihomo = self.moset.ret_ihomo()
            if imo%2==0:
                imo2 = ihomo - imo/2
            else:
                imo2 = ihomo + imo/2 + 1
            try:
                ene, occ = self.moset.ret_eo(imo2)
            except:
                print "\n ERROR: imo2 = %i"%imo2
                raise
            return "%s %.3f / %.3f %s"%(pref,ene,occ,postf)            

class mocoll_occ(mocoll):
    """
    Specify occupation threshold.
    """
    def __init__(self, occmin=0.01, occmax=2.0, mldfile=""):
        if mldfile == "":
            raise error_handler.MsgError("mldfile has to be specified for occupation screening!")
        
        self.mldfile = mldfile
        self.moset = lib_mo.MO_set_molden(mldfile)
        self.moset.read(lvprt=0)
        
        self.molist = []
        for imo, occ in enumerate(self.moset.occs):
            if occmin <= occ <= occmax:
                self.molist.append(imo)        

class jmol_options(input_options.write_options):
    def jmol_input(self):
        self.read_float('Cutoff value', 'cutoff', 0.05)
        
        #print ""
        self.choose_list(
            'Specification of the orbital indices to be plotted',
            'spec',
        [
            ('sten', 'Start and end indices'),
            ('frontier', 'Number of frontier orbitals'),
            ('occ', 'Occupation threshold')
        ]
        )
        #self.read_yn('Specification in terms of frontier orbitals', 'fr_mos')
        
        if self['spec'] == 'sten':
            self.read_int('First orbital index to be plotted', 'st_ind', 1)
            self.read_int('Last orbital index to be plotted',  'en_ind', 10)
        elif self['spec'] == 'frontier':
            self.read_int('Number of frontier orbitals',  'en_ind', 3)
        elif self['spec'] == 'occ':
            self.read_float('Minimal occupancy', 'occmin', 0.01)
            self.read_float('Maximal occupancy', 'occmax', 1.99)
        else:
            raise error_handler.ElseError(self['spec'], 'spec')
        
        self.read_yn('Use "rotate best" command (only available in Jmol 14)', 'rot_best')
        self.read_yn('Additional custom rotation of the molecule?', 'rot_custom')
        if self['rot_custom']:
            self.read_float('Rotation around the x-axis', 'rot_x', 0.)
            self.read_float('Rotation around the y-axis', 'rot_y', 0.)
            self.read_float('Rotation around the z-axis', 'rot_z', 0.)
        
        self.read_int('Width of images in output html file', 'width', 400)
        
def run():    
    print 'jmol_MOs.py [<mldfile> [<mldfile2> ...]]\n'
    
    mldfiles = sys.argv[1:]
    
    if len(mldfiles) == 0:
        print "No file specified, generating generic script"
        mldfiles = ['']
        pref = ''
    elif len(mldfiles) == 1:
        print "Analyzing the file:", mldfiles[0]
        pref = mldfiles[0] + '.'
    else:
        print "Analyzing the files:", mldfiles
        pref = 'multi.'
        
    jopt = jmol_options('jmol.in')
    jopt.jmol_input()
    
    jo = lib_file.wfile('%sjmol_orbitals.spt'%pref)
    ho = lib_file.htmlfile('%sorbitals.html'%pref)
    lo = lib_file.latexfile('%sorbitals.tex'%pref)
    
    ho.pre('Orbitals')
    lo.pre(None, graphicx=True)
    
    for mldfile in mldfiles:
        print 'Analyzing %s ...\n'%mldfile
        if jopt['spec'] == 'sten':
            moc = mocoll(jopt['st_ind'], jopt['en_ind'], mldfile)
        elif jopt['spec'] == 'frontier':
            moc = mocollf(jopt['en_ind'], mldfile)
        elif jopt['spec'] == 'occ':
            moc = mocoll_occ(jopt['occmin'], jopt['occmax'], mldfile)
        else:
            raise error_handler.ElseError(self['spec'], 'spec')
        
        moout = mo_output_jmol(moc, jopt)
        moout.output(jo)
        
        moh = mo_output_html(moc, jopt)
        moh.output(ho)
        
        mol = mo_output_tex(moc, jopt)
        mol.output(lo)
            
    jo.post(lvprt=1)
    if mldfiles == [""]:
        print "  -> Open the Molden-file in jmol and execute the commands contained in this file."
    else:
        print "  -> Now simply run \"jmol -n %s\" to plot all the orbitals.\n"%jo.name
    ho.post(lvprt=1)
    print "  -> View in browser."
    lo.post(lvprt=1)
    print "  -> Compile with pdflatex (or adjust first)."

if __name__=='__main__':
    import sys

    theo_header.print_header('Orbital plotting in Jmol')
    run()