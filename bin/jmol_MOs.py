#!/usr/bin/python
import lib_mo

class mo_output:
    """
    MO output in standard jmol format.
    """
    def __init__(self, moc):
        self.moc = moc
        
    def output(self):
        self.pre()
        self.print_mos()
        self.post()
        
    def pre(self, cutoff=0.05):
        if self.moc.mldfile != "":
            print '\nload "' + self.moc.mldfile + '" FILTER "nosort"'
        print "mo titleformat ''"
        print "background white\nmo fill\nmo cutoff %.3f\n"%cutoff
        
    def print_mos(self):
        for imo in self.moc.molist:
            print "mo %s"%self.moc.moname(imo)
            print "write image png \"%s\""%(self.moc.mopath(imo))
            
    def mopath(self, mo):
        return "MO_%s.png"%mo
            
    def post(self):
        pass
        
class mo_output_html(mo_output):
    """
    HTML file for visualizing the MOs created with jmol.
    """
    def __init__(self, moc, width=400):
        mo_output.__init__(self, moc)
        self.outfile = open("%sorbitals.html"%self.moc.ret_label(), "w")
        self.width = width
        
    def pre(self):
        print "\nWriting html file %sorbitals.html ..."%self.moc.ret_label(),
        
        self.outfile.write("<html>\n<head>\n<title>")
        self.outfile.write("Orbitals - %s"%self.moc.mldfile)
        self.outfile.write("</title>\n</head>\n<body>\n")
        self.outfile.write("<table>\n")
        
    def print_mos(self):
        for imo in self.moc.molist:
            if imo%2==0: self.outfile.write("<tr>\n")
            self.outfile.write("<td>")
            self.outfile.write("<img src=\"%s\""%(self.moc.mopath(imo)))
            self.outfile.write("border=\"1\" width=\"%i\">"%self.width)
            self.outfile.write(self.moc.mo_extra(imo,pref="<br> MO %i:"%imo))
            self.outfile.write("</td>\n")
            if imo%2==1: self.outfile.write("</tr>\n")
        if imo%2==0: self.outfile.write("</tr>\n")
        
    def post(self):
        self.outfile.write("</table>\n")
        self.outfile.write("</body>\n</html>\n")
        self.outfile.close()
    
        print "finished."
        
class mo_output_tex(mo_output):
    """
    tex file for visualizing the MOs created with jmol.
    """
    def __init__(self, moc, width=6., trim=[0.,0.,0.,0.]):
        mo_output.__init__(self, moc)
        self.outfile = open("%sorbitals.tex"%self.moc.ret_label(), "w")
        
        # take this out?
        self.igraphstr="["
        if not trim == [0.,0.,0.,0.]:
            self.igraphstr += "trim = %.2fcm %.2fcm %.2fcm %.2fcm, clip=true,"%(trim[0],trim[1], trim[2], trim[3])
        self.igraphstr += "width=%.2f cm]"%width

    def pre(self):
        print "\nWriting tex file %sorbitals.tex ..."%self.moc.ret_label(),
        
        self.outfile.write("\\documentclass[a4paper]{article}\n")
        self.outfile.write("\\usepackage[cm]{fullpage}\n\\usepackage{graphicx}\n\n")
        self.outfile.write("% trim: left, bottom, right, top\n")
        self.outfile.write("\\newcommand{\\incMO}{\\includegraphics[trim = 1.00cm 1.00cm 1.00cm 1.00cm, clip=true,width=6.00 cm]}\n\n")
        self.outfile.write("\\begin{document}\n\n")
        self.outfile.write("\\begin{figure}\n")
        self.outfile.write("\\caption{%s}\n"%self.moc.mldfile)
        self.outfile.write("\\begin{tabular}{c | c}\n")
        

    def print_mos(self):
        for imo in self.moc.molist:
            self.outfile.write("\\incMO{%s}"%(self.moc.mopath(imo)))
            # self.outfile.write("\\includegraphics%s{%s}"%(self.igraphstr, self.moc.mopath(imo)))
            if imo%2==0: self.outfile.write("&\n")
            else:
                self.outfile.write("\\\\\n")
                self.outfile.write(self.moc.mo_extra(imo-1,postf="&\n"))
                self.outfile.write(self.moc.mo_extra(imo,postf="\\\\\n"))

    def post(self):
        self.outfile.write("\\end{tabular}\n\\end{figure}\n\n")
        self.outfile.write("\\end{document}\n")
        self.outfile.close()
        
        print "finished."

        
class mocoll:
    def __init__(self, st_ind, en_ind, mldfile=""):
        self.molist = range(st_ind, en_ind+1)
        self.mldfile = mldfile
        
        if not mldfile=="":
            self.moset = lib_mo.MO_set_molden(mldfile)
            self.moset.read()
        
    def moname(self, imo):
        return str(imo)
        
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
            return self.mldfile.split('.')[0] + '_'
        
class mocollf(mocoll):
    """
    For frontier MOs.
    """
    def __init__(self, en_ind, mldfile=""):
        mocoll.__init__(self, 0, 2*en_ind-1, mldfile)
        
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

            


if __name__=='__main__':
    import sys

    ojmol = True
    ohtml = True
    otex  = True

    print "%s <st_ind> <en_ind> [<mldfile>]"%sys.argv[0]
    print "or:"
    print "%s -f <num_mo> [<mldfile>]\n"%sys.argv[0]

    if len(sys.argv)<3: sys.exit()

    if sys.argv[1] == '-f':
      fr_mos = True
    else:
      fr_mos = False
      st_ind = int(sys.argv[1])
    en_ind = int(sys.argv[2])

    if len(sys.argv)>=4: mldfiles = sys.argv[3:]
    else: mldfiles = [""]

    if ojmol:
        print " *** Copy into the jmol console:\n"
        for mldfile in mldfiles:
            if not fr_mos:
              moc = mocoll(st_ind, en_ind, mldfile)
            else:
              moc = mocollf(en_ind, mldfile)           
            moout = mo_output(moc)
            moout.output()
        print " \n *** jmol input finished"

    for mldfile in mldfiles:
        if not fr_mos:
          moc = mocoll(st_ind, en_ind, mldfile)
        else:
          moc = mocollf(en_ind, mldfile)                    
        if ohtml:
            moh = mo_output_html(moc)
            moh.output()
        if otex:
            mot = mo_output_tex(moc,width=3., trim=[5.0,1.0,5.0,3.0])
            mot.output()
            

