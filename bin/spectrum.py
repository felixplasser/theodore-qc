#!/usr/bin/python
"""
Create a convoluted spectrum from the oscillator strengths.
"""

import math
import theo_header, units, lib_file, input_options

do_plots = True
try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except:
    print "pylab/matplotlib not installed - plotting not possible"
    do_plots = False

class spec_options(input_options.write_options):
    def spec_input(self):
        self.read_int("Number of points in the spectrum", "npts", 200)
        self.read_float("Minimum energy in the plot (eV)", "emin", 2.0)
        self.read_float("Maximum energy in the plot (eV)", "emax", 8.0)
        self.read_float("FWHM broadening (eV)", "fwhm", 0.5)
        self.read_int("Lineshape: 1 - Lorentzian, 2 - Gaussian", "lineshape", 1)
        
        self.spec = spectrum(**self.opt_dict)
        
        self.read_str("Name of the file to analyze", "ana_file", "tden_summ.txt", True)
        
    def make_spec(self):
        sfile = lib_file.summ_file(self['ana_file'])
        header = sfile.ret_header()
        ddict  = sfile.ret_ddict()
        state_labels = sfile.ret_state_labels()
        
        for state in state_labels:
            try:
                f = ddict[state]['f']
            except KeyError:
                f = 0.           
            
            self.spec.add(f, ddict[state]['dE(eV)'])

        self.spec.normalize()
        
        self.spec.ascii_file()

        if do_plots:
            self.spec.plot(xunit='eV', pname='spectrum_eV.png')
            self.spec.plot(xunit='nm', pname='spectrum_nm.png')
            self.spec.plot(xunit='rcm', pname='spectrum_rcm.png')
        
# Code adapted from SHARC
class gauss:
    def __init__(self,fwhm):
        self.c=-4.*math.log(2.)/fwhm**2
        
    def ev(self,A,x0,x):
        return A*math.exp( self.c*(x-x0)**2)

class lorentz:
    def __init__(self,fwhm):
        self.c=0.25*fwhm**2
        
    def ev(self,A,x0,x):
        return A/( (x-x0)**2/self.c+1)

class spectrum:
    def __init__(self,npts,emin,emax,fwhm,lineshape):
      self.npts=npts
      (self.emin, self.emax) = (emin, emax)
      if lineshape==1:
          self.f=gauss(fwhm)
      elif lineshape==2:
          self.f=lorentz(fwhm)
      
      self.en = [emin + float(i)/self.npts*(emax-emin) for i in range(self.npts+1)]       # the energy grid needs to be calculated only once
      lamev = units.energy['nm'] * units.energy['eV']
      self.lam = [lamev / en for en in self.en]
      self.spec=[ 0. for i in range(self.npts+1) ]
      
      self.sticks=[] # list of pairs (A,x0), unit for x0: eV
      
    def add(self,A,x0):
        if A==0.:
            return
        
        self.sticks += [(A,x0)]
        
        for i in range(self.npts+1):
            self.spec[i]+=self.f.ev(A,x0,self.en[i])
    
    def normalize(self):
        smax = max(self.spec)
        print 'Normalizing the spectrum ...'
        print 'Maximum: % .5f'%smax
        for i, sp in enumerate(self.spec):
            self.spec[i] = sp / smax

        for i, stick in enumerate(self.sticks):
            self.sticks[i] = (stick[0] / smax, stick[1])

    def ascii_file(self, fname='spectrum.dat'):
        wf = lib_file.wfile(fname)
        
        wt = lib_file.asciitable(ncol=3)
        for i, en in enumerate(self.en):
            wt.add_row([en, self.spec[i], self.lam[i]])
        
        wf.write(wt.ret_table())
        wf.post(lvprt=1)
        
    def plot(self, xunit=1, pname='spectrum.png', lvprt=1):
        pylab.figure(figsize=(8,6))
                
        if xunit.lower() == 'ev':
            xlist = self.en
            pylab.xlabel('Energy (eV)')
            plot_sticks = self.sticks
            (xmin, xmax) = (self.emin, self.emax)
        elif xunit.lower() == 'nm':
            xfac = units.energy['nm'] * units.energy['eV']            
            xlist = self.lam
            #pylab.xlabel(r'$\lambda$') not working ...
            pylab.xlabel('Wavelength (nm)')
            plot_sticks = [(A, xfac/x0) for A,x0 in self.sticks]
            (xmin, xmax) = (xfac / self.emin, xfac / self.emax)
        elif xunit.lower() == 'rcm':
            xfac = 1. / units.energy['eV'] * units.energy['rcm']
            xlist = [en * xfac for en in self.en]
            pylab.xlabel('Wavenumber (1/cm)')            
            plot_sticks = [(A, x0 * xfac) for A,x0 in self.sticks]
            (xmin, xmax) = (self.emin * xfac, self.emax * xfac)
        else:
            raise error_handler.ElseError(xunit, 'xunit')
        pylab.ylabel('Absorption (norm.)')
        
        pylab.plot(xlist, self.spec, 'k-')
        
        for A,x0 in plot_sticks:
            pylab.plot([x0, x0], [-1., A], 'rx-')
        
        pylab.axis(xmin=xmin, xmax=xmax, ymin=0., ymax=1.)
        pylab.savefig(pname)
        
        if lvprt >= 1:
            print "Spectrum file %s created."%pname
    
if __name__ == '__main__':
    theo_header.print_header('Create a convoluted spectrum')
    
    sopt = spec_options('spectrum.in')
    sopt.spec_input()
    
    sopt.make_spec()