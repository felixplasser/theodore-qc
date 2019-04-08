"""
Parsing of files produced by different quantum chemical programs.
"""

from __future__ import print_function, division
from . import units, lib_mo, error_handler
import numpy
import os, struct


class file_parser_base:
    def __init__(self, ioptions):
        self.ioptions = ioptions

    def read(self):
        state_list = []

        return state_list

    def init_den(self, mos, rect=False):
        """
        Initialize an empty (transition) density matrix.
        rect=True specifies that only the the occ x (occ + virt) block is explicitly constructed.
           This allows to treat CIS like theories efficiently without many other changes in the code.
        """
        num_mo = mos.ret_num_mo()

        if not rect:
            return numpy.zeros([num_mo, num_mo])
        else:
            nocc  = mos.ret_ihomo() + 1
            nvirt = num_mo - nocc
            return numpy.zeros([nocc, num_mo])

    def parse_key(self, state, key, line, search_string, ind=-1, rfile=None, rmatrix=False, not_string=None):
        """
        Find search_string in the specified line and set it as state[key].
        If rfile is given, then the next line is searched for x,y,z components.
        If rmatrix==True, then a tensor is read from the next three lines.
        """
        translate_str = ',[]|' 
        if search_string in line:
            if not not_string is None:
                if not_string in line: return

            words = self.delete_chars(line.strip(search_string), translate_str)
            state[key] = float(words[ind])

            if not rfile is None:
                if rmatrix:
                    line = next(rfile)
                    try:
                        words=self.delete_chars(next(rfile), translate_str)
                        state['%sxx'%key] = float(words[-3])
                        state['%sxy'%key] = float(words[-2])
                        state['%sxz'%key] = float(words[-1])
                        words=self.delete_chars(next(rfile), translate_str)
                        state['%syx'%key] = float(words[-3])
                        state['%syy'%key] = float(words[-2])
                        state['%syz'%key] = float(words[-1])
                        words=self.delete_chars(next(rfile), translate_str)
                        state['%szx'%key] = float(words[-3])
                        state['%szy'%key] = float(words[-2])
                        state['%szz'%key] = float(words[-1])
                    except:
                        print(' Cannot parse matrix for %s of state %s'%(key, state['name']))
                else:
                    words=self.delete_chars(next(rfile), translate_str)
                    state['%sx'%key] = float(words[-3])
                    state['%sy'%key] = float(words[-2])
                    state['%sz'%key] = float(words[-1])

    def delete_chars(self, line, delete, lmap=None):
        """
        Delete the characters in <delete> from line, split, and optionally map.
        """
        retline = line
        for char in delete:
            retline = retline.replace(char, '')

        return retline.split()

    def sym_split(self, sym):
        """
        Split an MO "sym" label into the index and irrep.
        "51b1u" -> [51, "b1u"]
        """
        st = 1000
        if sym.find('a') != -1: st = min(st, sym.find('a'))
        if sym.find('b') != -1: st = min(st, sym.find('b'))
        if sym.find('e') != -1: st = min(st, sym.find('e'))
        if sym.find('t') != -1: st = min(st, sym.find('t'))

        if (st==1000): raise error_handler.MsgError("sym_split: %s"%sym)

        return (int(sym[:st]), sym[st:])

class configuration:
    def __init__(self):
        self.occ    = "1x" # label of the occupied (hole) orbital
        self.virt   = "1y" # label of the virtual (particle) orbital
        self.coeff  = 0.
        self.weight = 0.

    def read_ricc2_line(self, line):
        try:
            cols = line.split('|')

            words = cols[1].split()
            self.occ = words[0] + words[1]

            words = cols[2].split()
            self.virt = words[0] + words[1]

            words = cols[3].split()
            self.coeff  = float(words[0])
            self.weight = float(words[1])
        except:
            print("\n ERROR when parsing line:")
            print(line)
            raise

#--------------------------------------------------------------------------#
# Implementations
#--------------------------------------------------------------------------#

class file_parser_ricc2(file_parser_base):
    """
    Turbomole ricc2
    """
    def read(self, mos):
        state_list = self.ret_conf_ricc2(rfile=self.ioptions.get('rfile'))

        ihomo  = mos.ret_ihomo()

        for state in state_list:
            state['name'] = '%i(%s)%s'%(state['state_ind'], state['mult'], state['irrep'])
            state['tden'] = self.init_den(mos, rect=True)

            if self.ioptions.get('read_binary'):
                self.set_tden_bin(state, mos)
            else:
                self.set_tden_conf(state, mos)

        if self.ioptions.get('read_binary'):
               mos.symsort(self.ioptions['irrep_labels'], self.ioptions['Om_formula'])
               if self.ioptions['jmol_orbitals']:
                    print(" \nWARNING: jmol_orbitals not possible with read_binary. Use molden_orbitals instead!")
                    self.ioptions['jmol_orbitals'] = False

        return state_list

    def set_tden_conf(self, state, mos):
        """
        Set the transition density matrix elements acoording to the parsed configurations.
        """
        for conf in state['char']:
            # look for the orbital in the MO file. This makes sure that the procedure works even when the MOs are reordered because of symmetry.
            iocc  = mos.syms.index(conf.occ)

            #print conf.virt, mos.syms.index(conf.virt), mos.syms.index(conf.virt) - ihomo, conf.coeff
            ivirt = mos.syms.index(conf.virt)

            try:
                state['tden'][iocc,ivirt] = conf.coeff
            except:
                print("\nERROR setting %s->%s (%i->%i) transition"%(conf.occ, conf.virt, iocc+1, ivirt+1))
                print("Did you delete the line")
                print("       implicit core=   x virt=    x")
                print("  from the control file before running tm2molden?\n")
                raise

    def set_tden_bin(self, state, mos, lvprt=1):
        """
        Set the transition density matrix elements acoording to the binary files with state information.
        """
        iirrep=self.ioptions['irrep_labels'].index(state['irrep']) + 1
        smult=state['mult']
        istate=state['state_ind']
        CCfilen = 'CCRE0{0:->2d}{1:->3s}{2:->4d}'.format(iirrep, smult, istate)

        if lvprt >= 1: print('Reading binary file %s ...'%CCfilen)

        CCfile = open(CCfilen, 'rb')
        CCfile.read(8)
        method = struct.unpack('8s', CCfile.read(8))[0]
        CCfile.read(8)
        nentry = struct.unpack('i', CCfile.read(4))[0]
        CCfile.read(20)

        num_mo = mos.ret_num_mo()
        nocc  = mos.ret_ihomo() + 1
        nvirt = num_mo - nocc

        assert(nentry % nvirt == 0)
        nact = nentry // nvirt
        nfrzc = nocc - nact

        if lvprt >= 1:
            print('  Method: %s, number of entries: %i'%(method, nentry))
            print('  num_mo: %i, nocc: %i, nvirt: %i, nact: %i, nfrzc: %i'%(num_mo, nocc, nvirt, nact, nfrzc))
        if lvprt >= 2:
            print('  Syms:', mos.syms)

        if nfrzc > 0:
            raise error_handler.MsgError("""Frozen core orbitals detected in the Molden file!
In the case of read_binary=True, do not delete the line
       implicit core=   x virt=    x
from the control file.""")

        # write the collected data into the correct block of the 1TDM
        coeff = struct.unpack(nentry*'d', CCfile.read(nentry*8))
        state['tden'][nfrzc:nocc, nocc:num_mo] = numpy.reshape(coeff, [nact, nvirt])

        lbytes = struct.unpack('4s', CCfile.read(4))
        if lvprt >= 2:
            print('   Last four bytes:', lbytes)
        
        if not CCfile.read(1) == b'':
            raise error_handler.MsgError('parsing file %s'%CCfilen)

        if lvprt >= 3:
            print('parsed tden:')
            print(state['tden'])

    def ret_conf_ricc2(self, rfile='ricc2.out'):
        """
        Return information about configurations in a Turbomole calculation.
        """
        ret_list = []
        curr_osc = 0 # running index for oscillator strengths
        section  = '' # which section of the file
        lines = open(rfile,'r')
        while True: # loop over all lines
          try:
            line = next(lines)
          except StopIteration:
            print("Finished parsing file %s"%rfile)
            break
          else:
            if (' sym | multi | state' in line) and ('excitation energies' in line):
                # parse the multiplicity from here
                multis = []
                line = next(lines)
                line = next(lines)
                line = next(lines)
                while True:
                    line = next(lines)
                    if '---' in line: continue
                    if '===' in line: break

                    words = self.delete_chars(line, ['|'])

                    multis.append(words[1])

            elif 'Energy:   ' in line:
                ret_list.append({})
                words = line.split()
                ret_list[-1]['exc_en'] = float(words[3])

                for i in range(3): line = next(lines)

                words = line.split() # type: RE0
                ret_list[-1]['irrep'] = words[4]
                ret_list[-1]['state_ind'] = int(words[6])
                ret_list[-1]['mult'] = multis.pop(0)

                for i in range(3): line = next(lines)

                ret_list[-1]['char'] = []
                while True: # loop over information about this excitation
                    line = next(lines)

                    if '===' in line: break

                    conf = configuration()
                    conf.read_ricc2_line(line)
                    ret_list[-1]['char'].append(conf)

            elif 'oscillator strength (length gauge)' in line:
                words = line.split()
                if curr_osc < len(ret_list):
                    ret_list[curr_osc]['osc_str'] = float(words[5])
                elif curr_osc == len(ret_list):
                    print("\nWARNING: More oscillator strengths than singlet states!")
                    print("   This is ok if you are computing state-to-state properties.\n")
                curr_osc+=1

            elif curr_osc < len(ret_list):
                state = ret_list[curr_osc]
                self.parse_key(state, 'Tmux', line, 'xdiplen', 1)
                self.parse_key(state, 'Tmuy', line, 'ydiplen', 1)
                self.parse_key(state, 'Tmuz', line, 'zdiplen', 1)

            elif 'real representations' in line:
                self.ioptions['irrep_labels'] = line.split()[6:]

            elif 'COSMO-ADC(2) energy differences' in line:
                print('\nReading COSMO-ADC(2) energies ...')
                for i in range(4): line = next(lines)

                istate = 0
                while True: # loop over information about this excitation
                    line = next(lines)

                    if '---' in line: continue
                    if '===' in line: break

                    words = line.replace('|','').split()
                    irrep, mult, state_ind, exc_en = words[0], words[1], int(words[2]), float(words[-1])
                    state = ret_list[istate]
                    if state['irrep'] == irrep and state['mult'] == mult and state['state_ind'] == state_ind:
                        state['exc_en'] = exc_en
                    else:
                        print(('Data not matching:', state['irrep'], irrep, state['mult'], mult, state['state_ind'], state_ind))
                        raise error_handler.MsgError('COSMO-ADC: data not matching')

                    istate += 1

        return ret_list

#---

class file_parser_escf(file_parser_base):
    """
    Turbomole TDDFT
    """
    def read(self, mos):
        state_list = self.ret_conf_tddft(rfile=self.ioptions.get('rfile'))

        nocc={}
        nvirt={}

        for state in state_list:
            state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
            state['tden'] = self.init_den(mos, rect=True)

            occmap  = []
            virtmap = []
            for iorb, sym in enumerate(mos.syms):
                if state['irrep'] in sym:
                    occ = mos.occs[iorb]
                    if abs(occ-2.) < 1.e-4:
                        occmap.append(iorb)
                    elif abs(occ) < 1.e-4:
                        virtmap.append(iorb)
                    else:
                        print(" Error: invalid occupation!", occ)
                        exit(5)

            nocc=len(occmap)
            nvirt=len(virtmap)

            print("\n Considering state: %s"%state['name'])
            print("  Number of occupied orbitals:", nocc)
            print("  Number of virtual orbitals:", nvirt)
            print("  Tensor space dimension:", nocc*nvirt)

            if os.path.exists('sing_%s'%state['irrep']):
              readf = 'sing_%s'%state['irrep']
              print('  Reading information of singlet calculation from file %s'%readf)
            elif os.path.exists('trip_%s'%state['irrep']):
              readf = 'trip_%s'%state['irrep']
              print('  Reading information of triplet calculation from file %s'%readf)
            elif os.path.exists('ciss_%s'%state['irrep']):
              readf = 'ciss_%s'%state['irrep']
              print('  Reading information of TDA singlet calculation from file %s'%readf)
            elif os.path.exists('cist_%s'%state['irrep']):
              readf = 'cist_%s'%state['irrep']
              print('  Reading information of TDA triplet calculation from file %s'%readf)
            else:
              print('No file with information about the excited state (sing_a, trip_a, ...) found!')
              exit(7)

            # the file is parsed once for every state
            #   this is not efficient but also not critical
            curr_state = 0
            for line in open(readf):
                if 'tensor space dimension' in line:
                    words = line.split()
                    space_dim = int(words[-1])
                    assert(nocc*nvirt == space_dim)
                elif 'eigenvalue' in line:
                    words = line.split()
                    curr_state = float(words[0])
                    iiocc  = 0
                    iivirt = 0
                elif curr_state == state['state_ind']:
                    words = [line[0+20*i:20+20*i] for i in range(4)]
                    for word in words:
                        state['tden'][occmap[iiocc],virtmap[iivirt]] = float(word.replace('D','E'))

                        iivirt+=1
                        if iivirt == nvirt:
                            iivirt = 0
                            iiocc += 1
                            if iiocc == nocc:
                                curr_state = 0
                                break
                elif curr_state > state['state_ind']: break

        return state_list

    def ret_conf_tddft(self, rfile):
        rlines = open(rfile, 'r').readlines()[100:]
        ret_list = []
        occ_orb = False # section of the file
        for nr,line in enumerate(rlines):
            if 'excitation' in line and not 'vector' in line:
                ret_list.append({})
                words = line.split()
                ret_list[-1]['state_ind'] = int(words[0])
                #ret_list[-1]['irrep'] = file_handler.line_to_words(line)[2]
                ret_list[-1]['irrep'] = words[2]
                ret_list[-1]['tot_en'] = eval(rlines[nr+3][40:])
                ret_list[-1]['exc_en'] = eval(rlines[nr+7][40:])
                ret_list[-1]['osc_str'] = eval(rlines[nr+16][40:])
                ret_list[-1]['char'] = []
            elif 'occ. orbital' in line:
                occ_orb = True
                #print 'occ. orbital in line'
            elif occ_orb:
                words = line.split()
                #print words
                if len(words) > 0:
                    #print 'words > 0'
                    ret_list[-1]['char'].append({})
                    ret_list[-1]['char'][-1]['occ'] = words[0]+words[1]
                    ret_list[-1]['char'][-1]['virt'] = words[3]+words[4]
                    ret_list[-1]['char'][-1]['weight'] = float(words[-1])/100.
                    # escf does not print the phase of the coefficient
                else:
                    occ_orb = False

        return ret_list

#---

class file_parser_libwfa(file_parser_base):
    def __init__(self, ioptions):
        self.ioptions = ioptions
        self.excD = False
        self.excT = False

    def read(self):
        """
        Read the .om files as created with libwfa.
        These already contain OmAt
        """
        state_list = []

        basedir='.'
        #ist = 1
        for omfile in sorted(os.listdir(basedir)):
            suff = omfile.split('.')[-1]
            if suff!='om': continue

            (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile(os.path.join(basedir,omfile))
            if typ == None: continue

            state_list.append({})

            # Here several try statements are used to check for different versions of the output
            #    Not very elegant but it should work ...
            try:
                multl, iirrep, ist = typ.split('_')
            except:
                state_list[-1]['name'] = typ
            else:
                try:
                    state_list[-1]['irrep'] = self.ioptions.get('irrep_labels')[int(iirrep)]
                except:
                    state_list[-1]['irrep'] = iirrep

                state_list[-1]['mult'] = {'singlet':'(1)',
                                          'triplet':'(3)',
                                          'any':'(-)'}[multl]

                state_list[-1]['state_ind'] = int(ist)

                state_list[-1]['name'] = '%i%s%s'%(state_list[-1]['state_ind'],
                                                    state_list[-1]['mult'],
                                                    state_list[-1]['irrep'])

            #state_list[-1]['char'] = []
            state_list[-1]['exc_en'] = exctmp * units.energy['eV']
            state_list[-1]['osc_str'] = osc

            state_list[-1]['Om']   = om_at.sum()
            state_list[-1]['OmAt'] = om_at

          #ist += 1

        return state_list

    def rmatfile(self,fname):
        print("Reading: %s ..."%fname)

        try:
            rfile=open(fname,'r')
        except IOError:
            print("\n WARNING: could not open %s."%fname)
            print("Did the calculation converge?\n")
            return None, None, None, None, None, None

        line=rfile.readline()
        words=line.split()

        if '<-->' in line:
            # ccman2 output
            typ = words[2]
            excen = float(words[3]) if len(words) >= 4 else  0.
            osc   = float(words[4]) if len(words) >= 5 else -1.
        else:
            typ = words[0]
            excen = float(words[1]) if len(words) >= 2 else  0.
            osc   = float(words[2]) if len(words) >= 3 else -1.

        line=rfile.readline()
        words=line.split()
        dima=int(words[1])
        dimb=int(words[2])

        #print " Dimensions: %i x %i"%(dima, dimb)

        outarr = numpy.zeros([dima,dimb])

        ia = 0
        ib = 0
        while 1:
          line = rfile.readline()
          if not line:
            break
          vals = [float(word) for word in line.split()]

          for val in vals:
            outarr[ib, ia] = val
            ib += 1
            if ib == dimb:
              ia += 1
              ib = 0
              if ia == dima:
                print(" All values read in")

        return typ, excen, osc, dima, dimb, outarr

    def parse_line(self, state, line, rfile=None):
        if 'Exciton analysis of the difference density matrix' in line:
            self.excD = True
            self.excT = False
        elif 'Exciton analysis of the transition density matrix' in line:
            self.excD = False
            self.excT = True

        self.parse_keys(state, self.excD, self.excT, line, rfile)

    def parse_keys(self, state, exc_diff, exc_1TDM, line, rfile=None):
        self.parse_key(state, 'mu', line, 'Total dipole')
        self.parse_key(state, 'mux', line, 'Dip. moment [a.u.]', -3)
        self.parse_key(state, 'muy', line, 'Dip. moment [a.u.]', -2)
        self.parse_key(state, 'muz', line, 'Dip. moment [a.u.]', -1)
        self.parse_key(state, 'Tmux', line, 'Trans. dip. moment [a.u.]', -3)
        self.parse_key(state, 'Tmuy', line, 'Trans. dip. moment [a.u.]', -2)
        self.parse_key(state, 'Tmuz', line, 'Trans. dip. moment [a.u.]', -1)
        self.parse_key(state, 'r2x', line, '<r^2>', -3, not_string='Total')
        self.parse_key(state, 'r2y', line, '<r^2>', -2, not_string='Total')
        self.parse_key(state, 'r2z', line, '<r^2>', -1, not_string='Total')
        self.parse_key(state, 'r2', line, 'Total <r^2>')
        self.parse_key(state, 'nu', line, 'Number of unpaired electrons:', 2)
        self.parse_key(state, 'nunl', line, 'Number of unpaired electrons')
        self.parse_key(state, 'p', line, 'Number of detached / attached electrons')
        self.parse_key(state, 'Om_', line, 'Sum of SVs')
        self.parse_key(state, 'PRNTO', line, 'PR_NTO')
        self.parse_key(state, 'PRD', line, 'PR_D', 6)
        self.parse_key(state, 'PRA', line, 'PR_A')
        self.parse_key(state, '2P', line, 'Two-photon absorption cross-section [a.u.]', rfile=rfile, rmatrix=True)
        self.parse_key(state, 'S_HE', line, 'Entanglement entropy')
        self.parse_key(state, 'Z_HE', line, 'Nr of entangled states')
        self.parse_key(state, 'mu', line, 'Dipole moment [D]', rfile=rfile)
        self.parse_key(state, 'sigR', line, 'RMS size of the density', rfile=rfile)

        if exc_diff:
            self.parse_key(state, 'sigD', line, 'Hole size', rfile=rfile)
            self.parse_key(state, 'sigA', line, 'Electron size', rfile=rfile)
            self.parse_key(state, 'dD-A', line, '|<r_e - r_h>|')

        if exc_1TDM:
            self.parse_key(state, 'dexc', line, 'RMS electron-hole separation', rfile=rfile)
            self.parse_key(state, 'sigH', line, 'Hole size', rfile=rfile)
            self.parse_key(state, 'sigE', line, 'Electron size', rfile=rfile)
            self.parse_key(state, 'dH-E', line, '|<r_e - r_h>|')
            self.parse_key(state, 'COV', line, 'Covariance(r_h, r_e) [Ang^2]')
            self.parse_key(state, 'Corr', line, 'Correlation coefficient')

class file_parser_qcadc(file_parser_libwfa):
    """
    Parse information from qchem.out in addition to the .om file.
    """
    def read(self):
        state_list = []

        basedir='.'

        exc_diff = False
        exc_1TDM = False
        self.irrep_labels = None
        rf = open(self.ioptions.get('rfile'), 'r')
        for line in rf:
            words = line.split()
            if 'Irreducible representations in point group:' in line:
                self.irrep_labels = line.lstrip('Irreducible representations in point group:').split()

            elif ' Term symbol' in line:
                if self.irrep_labels == None:
                    print("\n   WARNING: irrep labels not found in %s"%self.ioptions.get('rfile'))
                    print("   Use 'adc_print 2' or enter them into %s"%self.ioptions.ifile)
                    print("   Using info from %s: irrep_labels="%self.ioptions.ifile, self.ioptions.get('irrep_labels'))
                    print()
                    self.irrep_labels = self.ioptions.get('irrep_labels')

                state_list.append({})

                state_list[-1]['state_ind'] = int(words[2])
                state_list[-1]['mult']      = words[3]
                state_list[-1]['irrep']     = words[4]
                state_list[-1]['name']      = '%s%s%s'%(words[2], words[3], words[4])

                om_filen = self.om_file_name(state_list[-1])
                (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile(os.path.join(basedir,om_filen))
                if typ == None:
                    continue

                state_list[-1]['exc_en'] = exctmp * units.energy['eV']
                state_list[-1]['osc_str'] = osc

                state_list[-1]['Om']   = om_at.sum()
                state_list[-1]['OmAt'] = om_at

            elif 'MP(2) Summary' in line:
                state_list.append({})
                state_list[-1]['name'] = 'gr-st'
                state_list[-1]['exc_en'] = 0.
                state_list[-1]['state_ind'] = 0

            elif ' Excitation energy:' in line:
                exc_chk = float(words[2])

                if not 'exc_en' in state_list[-1]:
                    state_list[-1]['exc_en'] = exc_chk

                if abs(exc_chk - state_list[-1]['exc_en']) > exc_chk * 1.e-4:
                    print(exc_chk, state_list[-1]['exc_en'])
                    raise error_handler.MsgError("Excitation energies do not match")

            elif 'Exciton analysis of the difference density matrix' in line:
                exc_1TDM = False
                if len(state_list) > 0: exc_diff = True

            elif 'Exciton analysis of the transition density matrix' in line:
                exc_diff = False
                if len(state_list) > 0: exc_1TDM = True

            elif 'Transition Summary' in line:
                break

            if len(state_list) > 0:
                self.parse_keys(state_list[-1], exc_diff, exc_1TDM, line, rf)

        rf.close()

        # Post-processing
        # Pre-factor for 2P cross-section
        pre2P = numpy.pi**3 * units.constants['c0']**(-2) * units.tpa['GM']
        print("pre2P, GM: ", pre2P, units.tpa['GM'])
        # This is numerically the same as [JCP, 146, 174102]:
        #   numpy.pi**3 / u.c0 * u.cm**5 / 2.9979E10 * 10**50
        for state in state_list:
            # 2P cross-section pre-factor in GM
            if '2P' in state:
                state['GM'] = pre2P * state['2P']/30 * (state['exc_en'] / units.energy['eV'])
                state['2P'] *= .000001
            if 'Tmux' in state:
                state['Tmu'] = (state['Tmux']**2 + state['Tmuy']**2 + state['Tmuz']**2)**.5

        return state_list

    def om_file_name(self, state):
        """
        Construct the name of the .om file
        This routine parses different versions of Q-Chem and libwfa.
        """
        multo = {'(1)':'singlet',
                 '(3)':'triplet',
                 '(-)':'any'}[state['mult']]

        # subtract 1 for ground state irrep
        state_indo = state['state_ind']
        if state['mult'] == '(1)':
            if state['irrep'].lower() in ['a', 'ag', 'a1', "a'"]:
                state_indo -= 1

        state['fname'] = '%s_%s_%i'%(multo, state['irrep'], state_indo)
        omname = '%s_ctnum_atomic.om'%state['fname']

        if not os.path.exists(omname):
            print((' WARNING: did not find %s,' % omname), end=' ')
            irrepo = self.irrep_labels.index(state['irrep'])

            # subtract 1 for ground state irrep
            state_indo = state['state_ind'] - 1*(irrepo==0)*(state['mult'] == '(1)')

            state['fname'] = '%s_%i_%i'%(multo, irrepo, state_indo)
            omname = '%s_ctnum_atomic.om'%state['fname']

            print('using %s instead.'%omname)

        return omname

class file_parser_qctddft(file_parser_libwfa):
    def read(self, mos):
        """
        Read the X vector from standard output. Y is discarded.
        """
        state_list = []
        exc_diff = exc_1TDM = tdread = libwfa = False
        istate = 1

        if self.ioptions.get('TDA'):
            ststr = 'TDDFT/TDA Excitation Energies'
        else:
            ststr = 'TDDFT Excitation Energies'

        print("Parsing %s for %s ..."%(self.ioptions.get('rfile'), ststr))

#        for line in open(self.ioptions.get('rfile'), 'r'):
        rfileh = open(self.ioptions.get('rfile'),'r')
        while True: # loop over all lines
            try:
                line = next(rfileh)
            except StopIteration:
              print("Finished parsing file %s"%self.ioptions.get('rfile'))
              break

            if ststr in line:
                tdread = True
            elif 'TDDFT calculation will be performed' in line:
                tdread = False
            elif 'Timing summary' in line:
                tdread = False
            elif 'Excited State Analysis' in line:
                tdread = False
                libwfa = True
                if len(state_list) == 0:
                    errstr = "No excitation energy parsed!"
                    errstr+= "\n   Please, set 'TDA=True' if this was a TDA calculation."
                    raise error_handler
            elif 'SA-NTO Decomposition' in line:
                libwfa = False
            elif 'Welcome to Q-Chem' in line:
                if len(state_list) > 0:
                    print("\n WARNING: found second Q-Chem job!\n   Deleting everything parsed so far.\n")

                state_list = []
                exc_diff = exc_1TDM = tdread = libwfa = False
                istate = 1

            if tdread:
                words = line.replace(':','').split()
                if 'Excited state' in line:
                    state_list.append({})

                    state_list[-1]['state_num'] = int(words[2])
                    state_list[-1]['exc_en'] = float(words[-1])

                    line = next(rfileh)
                    line = next(rfileh)
                    words = line.split()

                    if words[0] == 'Multiplicity:':
                        state_list[-1]['mult'] = '(%s)'%words[1][0]
                    else:
                        state_list[-1]['mult'] = '(-)'

                    state_list[-1]['name'] = 'es_%i%s'%(state_list[-1]['state_num'], state_list[-1]['mult'])

                    if self.ioptions['read_libwfa']:
                        om_filen = 'es_%i_ctnum_atomic.om'%state_list[-1]['state_num']
                        (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile(om_filen)

                        if not typ==None:
                            state_list[-1]['Om']   = om_at.sum()
                            state_list[-1]['OmAt'] = om_at
                    else:
                        state_list[-1]['tden'] = self.init_den(mos, rect=True)

                elif 'Strength' in line:
                    state_list[-1]['osc_str'] = float(words[-1])

                elif 'amplitude' in line and not self.ioptions['read_libwfa']:
                    # ignore the Y vector.
                    #    Otherwise the Y would go into the virt-occ block!
                    if 'Y:' in line: continue

                    awords = self.delete_chars(line, ['X:', 'Y:', 'D', 'V', '(', ')', '-->'])

                    iocc = int(awords[0]) - 1
                    ivirt = int(awords[1]) + mos.ret_ihomo()

                    coeff =  float(awords[4])

                    state_list[-1]['tden'][iocc, ivirt] += coeff

            if libwfa:
                if 'Excited state' in line:
                    istate = int(line.split()[-1].replace(':',''))
                elif 'Exciton analysis of the difference density matrix' in line:
                    exc_1TDM = False
                    exc_diff = True

                elif 'Exciton analysis of the transition density matrix' in line:
                    exc_diff = False
                    exc_1TDM = True

                self.parse_keys(state_list[istate-1], exc_diff, exc_1TDM, line)

        rfileh.close()
        return state_list

#---

class file_parser_col(file_parser_base):
    def read_iwfmt(self, dens, filen, fac = 1.):
        """
        Read output from iwfmt.x for a 1-particle density file.
        """
        eref = eexc = 0.

        in_data=False
        #for line in open(filen,'r').readlines():
        rfile = open(filen, 'r')
        while True:
            try:
                line=next(rfile)
            except StopIteration:
                break
            if ' 0.000000000000E+00' in line:
                if len(line.split()) == 1:
                    words = last_line.split()
                    (num,lab1,ibvtyp,itypea,itypeb,ifmt,last,nipv) = [int(word) for word in words]
                    if itypea==0 and itypeb==7:
                        print(' Reading symmetric density ...')
                        in_data=True
                        sym = 1
                    elif itypea==2 and itypeb==9:
                        print(' Reading antisymmetric density ...')
                        in_data=True
                        sym = -1
                    else:
                        print('Warning: this section does not contain a density')
                        print('itypea: %i, itypeb: %i'%(itypea,itypeb))
                        print(last_line)
                        in_data=False
            elif '-1025   -1025      -1' in line:
                line=next(rfile)
                (eexc, eref, ecore) = [float(word) for word in line.split()]
            elif '    ' in line:
                pass
            elif in_data:
                words = line.split()

                val = float(words[0])
                i = int(words[1])-1
                j = int(words[2])-1

                if 0.2 < abs(val) < 1.9999: print("(i,j)=(%2i,%2i), val=%6.3f"%(i,j,val))
                dens[j,i] += fac * val

                if i != j: dens[i,j] += fac * sym * val

            last_line = line

        return eref, eexc

    def read_sifs(self):
        """
        Directly read a SIFS file.
        This would require a small Fortran program with the required routines.
            -> Should be compiled with Columbus.
        Interface with f2py
        """
        raise error_handler.NIError()

    def read_trncils(self, state, mos, filen):
        """
        Read output from transci.x for a 1-particle density file.
        """
        # change the format of the MO labels
        mos.syms2 = [sym.lower().replace('_','') for sym in mos.syms]
        if not self.ioptions['ncore'] == {}:
            for isym, sym2 in enumerate(mos.syms2):
                imo, irr = self.sym_split(sym2)
                new = "%i%s"%(imo-self.ioptions['ncore'][irr], irr)
                mos.syms2[isym] = new

        tmp = filen.replace('LISTINGS/','').replace('state','').replace('drt','').replace('trncils.FROM','').split('TO')
        ir_st_ref = tmp[0].split('.')
        ir_st_exc = tmp[1].split('.')

        state['irrep'] = self.ioptions.get('irrep_labels')[int(ir_st_exc[0]) - 1]
        state['state_ind'] = int(ir_st_exc[1])
        state['name'] = '%s.%i-%i'%(state['irrep'], int(ir_st_ref[1]), state['state_ind'])

        if ir_st_ref==ir_st_exc:
            state['sden'] = self.init_den(mos)
            dfac = 1.
            mat = state['sden']
        else:
            state['tden'] = self.init_den(mos)
            dfac = 1. / numpy.sqrt(2)
            mat = state['tden']

        rfile = open(filen,'r')
        while True:
            try:
                line=next(rfile)
            except StopIteration:
                break
            if ' symm-trans density matrix:' in line:
                print('file %s, reading:'%filen)
                print(line.strip('\n'))
                line=next(rfile)
                self.read_block_mat(mat, mos, rfile, sym=1, dfac=dfac)
            elif 'asymm-trans density matrix:' in line:
                print('reading ...')
                print(line[:-1])
                line=next(rfile)
                self.read_block_mat(mat, mos, rfile, sym=-1, dfac=dfac)
            elif 'Transition energy:  ' in line:
                words=line.split() # use the line with "a.u." because it is given with more digits
                state['exc_en']=float(words[2]) * units.energy['eV']
                next(rfile)
            elif 'Oscillator strength : ' in line:
                words = line.split()
                state['osc_str']=float(words[-1])
        rfile.close()

    def read_block_mat(self, mat, mos, rfile, sym, dfac=1.):
        """
        Parse the block matrix output in the listing file.
        This is not the cleanest routine but it seems to work ...
        """
        while(1):
            words=next(rfile).replace('MO','').split()

            if 'density' in words: # skip two lines in case the next symmetry block is coming
                line=next(rfile)
                words=next(rfile).replace('MO','').split()

            try:
                head_inds = [mos.syms2.index(sym2) for sym2 in words]
            except:
                print(" ERROR reading:", words)
                print("syms2: ", mos.syms2)
                print()
                raise

            words=next(rfile).replace('MO','').split()
            while(len(words)>0): # loop over a block with constant header labels
                if 'integral' in words: break
                try:
                    left_ind = mos.syms2.index(words[0])
                except:
                    print('\n ERROR parsing: ')
                    print(words)
                    raise

                for i, word in enumerate(words[1:]):
                    val = float(word)
                    if val*val > 0.01: print("(%2i->%2i)-(%s->%s), val=% 8.4f"%\
                       (head_inds[i],left_ind,mos.syms2[head_inds[i]],mos.syms2[left_ind],val))

                    mat[head_inds[i],left_ind] += dfac * val
                    if left_ind != head_inds[i]:
                        mat[left_ind,head_inds[i]] += dfac * sym*val

                words=next(rfile).replace('MO','').split()

            if 'integral' in words: break

class file_parser_col_mrci(file_parser_col):
    def read(self, mos):
        if self.ioptions['s_or_t'] == 's':
            raise error_handler.MsgError('analyze_sden.py not implemented for col_mrci! Use "nos" instead.')
            # it is actually possible to have bra=ket transitions in transmomin ...

        state_list = []

        for lfile in sorted(os.listdir('LISTINGS')):
            if not 'trncils' in lfile: continue

            print("Reading %s ..."%lfile)
            state = {}
            self.read_trncils(state, mos, 'LISTINGS/%s'%lfile)
            if 'tden' in state:
                state_list.append(state)

        return state_list

    def read_fcd(self, mos):
        state_list = [{}, {}, {}]

        ist1 = self.ioptions['state_pair'][0]
        ist2 = self.ioptions['state_pair'][1]

        fn1  = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist1, ist1)
        fn2  = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist2, ist2)
        fn12 = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist1, ist2)
        fn21 = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist2, ist1)

        self.read_trncils(state_list[0], mos, fn1)
        self.read_trncils(state_list[1], mos, fn2)
        if os.path.exists(fn12):
            self.read_trncils(state_list[2], mos, fn12)
        elif os.path.exists(fn21):
            self.read_trncils(state_list[2], mos, fn21)
        else:
            print("Not using transition charges.")
            state_list[2] == None

        return state_list

class file_parser_col_mcscf(file_parser_col):
    def read(self, mos):
        state_list = []

        for lfile in sorted(os.listdir('WORK')):
            # Find the suitable files. This could also be done with regexps ...
            if not '.iwfmt' in lfile: continue
            if not 'mcsd1fl' in lfile: continue

            if self.ioptions['s_or_t'] == 't':
                if not '-' in lfile: continue

                print("Reading %s ..."%lfile)
                state_list.append({})
                self.read_mc_tden(state_list[-1], mos, lfile)

            elif self.ioptions['s_or_t'] == 's':
                if '-' in lfile: continue

                print("Reading %s ..."%lfile)
                state_list.append({})
                self.read_mc_sden(state_list[-1], mos, lfile)

        if len(state_list) == 0:
            raise error_handler.MsgError('No density file found! Did you run write_den.bash?')

        return state_list

    def read_fcd(self, mos):
        state_list = [{}, {}, {}]

        ist1 = self.ioptions['state_pair'][0]
        ist2 = self.ioptions['state_pair'][1]

        fn1  = "mcsd1fl.drt1.st%2.2i.iwfmt"%ist1
        fn2  = "mcsd1fl.drt1.st%2.2i.iwfmt"%ist2
        fn12 = "mcsd1fl.drt1.st%2.2i-st%2.2i.iwfmt"%(ist1, ist2)
        fn21 = "mcsd1fl.drt1.st%2.2i-st%2.2i.iwfmt"%(ist2, ist1)

        self.read_mc_sden(state_list[0], mos, fn1)
        self.read_mc_sden(state_list[1], mos, fn2)
        try:
            self.read_mc_tden(state_list[2], mos, fn12)
        except IOError:
            print("File %s not found, trying %s..."%(fn12, fn21))
            try:
                self.read_mc_tden(state_list[2], mos, fn21)
            except IOError:
                print("Not using transition charges.")
                state_list[2] == None

        return state_list

    def read_mc_tden(self, state, mos, filen):
        tmp = filen.replace('mcsd1fl.drt','').replace('st','').replace('.iwfmt','').split('.')
        state['irrep'] = self.ioptions.get('irrep_labels')[int(tmp[0]) - 1]
        state['name'] = '%s.%s'%(state['irrep'], tmp[1])
        state['tden'] = self.init_den(mos)

        # read symmetric part
        (eref, eexc) = self.read_iwfmt(state['tden'], 'WORK/'+filen, fac=1/numpy.sqrt(2.))
        state['exc_en'] = (eexc - eref) * units.energy['eV']

        # read antisymmetric part
        self.read_iwfmt(state['tden'], 'WORK/'+filen.replace('mcsd1fl', 'mcad1fl'), fac=1/numpy.sqrt(2.))

    def read_mc_sden(self, state, mos, filen):
        tmp = filen.replace('mcsd1fl.drt','').replace('st','').replace('.iwfmt','').split('.')
        state['irrep'] = self.ioptions.get('irrep_labels')[int(tmp[0]) - 1]
        state['state_ind'] = int(tmp[1])
        state['name'] = '%s.%i'%(state['irrep'], state['state_ind'])
        state['sden'] = self.init_den(mos)

        (eref, eexc) = self.read_iwfmt(state['sden'], 'WORK/'+filen)
        state['exc_en'] = eexc

class file_parser_terachem(file_parser_base):
    """
    Read terachem TDDFT job.
    """
    def read(self, mos):
        state_list = []

        rfileh = open(self.ioptions['rfile'], 'r')
        while True: # loop over all lines
            try:
                line = next(rfileh)
            except StopIteration:
              print("Reached end of file %s"%self.ioptions.get('rfile'))
              break

            if 'Largest CI coefficients' in line:
                state_list.append({})
                state = state_list[-1]
                state['tden'] = self.init_den(mos, rect=True)

                while True:
                    line = next(rfileh)
                    words = line.split()
                    if len(words) == 0: break

                    iocc  = int(words[0]) - 1
                    avirt = int(words[2]) - 1
                    if ('X' in line) and ('Y' in line): # rpa case
                        coeff = float(line.split(':')[2].split('Y')[0])
                    else:
                        coeff = float(words[-1])

                    state['tden'][iocc, avirt] = coeff

            elif 'Final Excited State Results' in line:
                line = next(rfileh)
                line = next(rfileh)
                line = next(rfileh)
                for state in state_list:
                    line = next(rfileh)
                    words = line.split()

                    state['state_ind'] = int(words[0])
                    state['name'] = 'A' + words[0]
                    state['exc_en'] = float(words[2])
                    state['osc_str'] = float(words[3])

        rfileh.close()

        return state_list

class file_parser_tddftb(file_parser_base):
    """
    Read DFTB+ TDDFTB job.
    Author: Ljiljana Stojanovic
    """
    def read(self, mos):
        state_list = []  # check if state_list is already allocated
        occ = []
        virt = []
        words = []
        Coeff = False
        # read the order of iocc, avirt pairs from the SPX.DAT file. It corresponds to the order in which X+Y coeffs are written in XplusY.DAT file
        spxfileh = open('SPX.DAT', 'r') # spxfile is SPX.DAT
        for line in spxfileh:
            words = line.split()
            if len(words) == 6:
               occ.append(int(words[3]))
               virt.append(int(words[5]))
        spxfileh.close()

        words = []
        rfileh = open('XplusY.DAT', 'r') # read X+Y coeffs from XplusY.DAT
        k = 0
        for line in rfileh:
            words = line.split()
            if 'S' in line:
                 k += 1
                 Coeff = True
                 xply = []
                 state_list.append({})
                 state = state_list[-1]
                 state['tden'] = self.init_den(mos, rect=True)
            elif Coeff:
               try:
                 if words[1] != 'S':
                     i = -1
                     for x in range (0,len(words)):
                         i += 1
                         xply.append(float(words[i]))
                 else: continue
                 j = -1
                 for x in range (0,len(xply)-1):
                     j += 1
                     iocc = occ[j]-1
                     avirt = virt[j]-1
                     state['tden'][iocc, avirt] = xply[j]
               except:
                 if len(words) == 3 and words[1]=='S': continue
               state['tden']=state['tden']/numpy.linalg.norm(state['tden'])

        rfileh.close()

        excfile = open(self.ioptions['rfile'], 'r') # rfile is EXC.DAT
        line = next(excfile)
        line = next(excfile)
        line = next(excfile)
        line = next(excfile)
        line = next(excfile)
        state_ind = 0
        for state in state_list:
            line = next(excfile)
            words = line.split()
            if len(words) == 8:
               state_ind += 1
               state['state_ind'] = state_ind
               state['name'] = str(state_ind) + 'A'
               state['exc_en'] = float(words[0])
               state['osc_str'] = float(words[1])

        excfile.close()

        return state_list

class file_parser_adf(file_parser_base):
    """
    Read ADF TDDFT using the TAPE21 file.
    """
    def read(self, mos):
        try:
            from scm.plams import KFFile
        except ImportError:
            from kf import kffile as KFFile

        rfile = KFFile(self.ioptions['rfile'])
        try:
            nmo = int(rfile.read('A','nmo_A'))
            nelec = int(rfile.read('General','electrons'))
        except:
            print(("\n  ERROR: reading TAPE21 file (%s)!\n"%self.ioptions['rfile']))
            raise
        assert nelec%2==0, "Odd number of electrons not supported"

        group = rfile.read('Geometry','grouplabel')[0]
        if not group in ['NOSYM', 'C1', 'c1']:
            print(("grouplabel: " + group))
            raise error_handler.MsgError('No support for symmetry')

        nocc = nelec / 2
        assert nocc == mos.ret_ihomo() + 1
        nvirt = nmo - nocc
        assert nvirt == mos.ret_num_mo() - nocc

        try:
            self.nsing = len(rfile.read('All excitations','All Sing-Sing excitations'))
            excss = rfile.read('Excitations SS A', 'excenergies') * units.energy['eV']
            oscss = rfile.read('Excitations SS A', 'oscillator strengths')
        except TypeError:
            self.nsing = 0
        try:
            self.ntrip = len(rfile.read('All excitations','All Sing-Trip excitations'))
            excst = rfile.read('Excitations ST A', 'excenergies') * units.energy['eV']
        except TypeError:
            self.ntrip = 0

        state_list = [{} for istate in range(self.nsing+self.ntrip)]
        istate = 0
        for ising in range(self.nsing):
            state = state_list[istate]

            state['state_ind'] = ising + 1
            state['exc_en'] = excss[ising]
            state['osc_str'] = oscss[ising]
            state['mult'] = 1
            state['name'] = '%i(%i)A'%(state['state_ind'], state['mult'])

            state['tden'] = self.init_den(mos, rect=True)
            eigen = rfile.read('Excitations SS A','eigenvector %s'%(istate+1))
            state['tden'][:,nocc:nmo] = eigen.reshape(nocc, nvirt)
            istate += 1

            # print-out
            print((state['name']))
            tden = state['tden']
            for i in range(len(tden)):
                for j in range(len(tden[0])):
                    val = tden[i, j]
                    if val*val > 0.1:
                        print(("(%i -> %i) % .4f"%(i+1,j+1,val)))

        for itrip in range(self.ntrip):
            state = state_list[istate]

            state['state_ind'] = itrip + 1
            state['exc_en'] = excst[itrip]
            state['mult'] = 3
            state['name'] = '%i(%i)A'%(state['state_ind'], state['mult'])

            state['tden'] = self.init_den(mos, rect=True)
            eigen = rfile.read('Excitations ST A','eigenvector %s'%(istate+1 - self.nsing))
            state['tden'][:,nocc:nmo] = eigen.reshape(nocc, nvirt)
            istate += 1

            # print-out
            print((state['name']))
            tden = state['tden']
            for i in range(len(tden)):
                for j in range(len(tden[0])):
                    val = tden[i, j]
                    if val*val > 0.1:
                        print(("(%i -> %i) % .4f"%(i+1,j+1,val)))

        rfile.close()

        return state_list

class file_parser_nos(file_parser_base):
    """
    Interpret the MO-file as a diagonal density.
    """
    def read(self, mos):
        state_list = []

        for no_file in self.ioptions['ana_files']:
            state_list.append({})
            self.read_no_file(state_list[-1], mos, no_file)

        for istate, state in enumerate(state_list):
            state['exc_en'] = float(istate + 1) # set fake excitation energy
            state['state_num'] = istate + 1

            state['fname'] = self.ioptions['ana_files'][istate]

            # perform some simplifications to the filename
            # for Columbus:
            tmp_name = state['fname'].replace('MOLDEN/','').replace('./','').replace('molden_no_','').rstrip('.sp')
            tmp_name = tmp_name.replace('state', 'S').replace('drt', 'D')
            # for Q-Chem
            tmp_name = tmp_name.replace('NOs/','').replace('.mo', '')

            state['name'] = tmp_name

        return state_list

    def read_ref_nos(self, state, ref_nos):
        """
        Read the reference MOs.
        ### This is currently not used! ###
        """
        state['sden'] = numpy.diag(ref_nos.occs)

        if self.ioptions['unpaired_ana']:
            nu_list = [min(occ, 2.-occ) for occ in ref_nos.occs]
            state['nu'] = sum(nu_list)
            state['nu_den'] = numpy.diag(nu_list)

            nunl_list = [occ*occ*(2-occ)*(2-occ) for occ in ref_nos.occs]
            state['nunl'] = sum(nunl_list)
            state['nunl_den'] = numpy.diag(nunl_list)

    def read_no_file(self, state, ref_mos, no_file):
        """
        Read information from a secondary NO file.
        """
        nos = lib_mo.MO_set_molden(file=no_file)
        nos.read()
        nos.compute_inverse()
        if self.ioptions['rd_ene']:
            nos.set_ens_occs()

        if not self.ioptions['ignore_irreps'] == []:
            print("  Ignoring irreps: ", self.ioptions['ignore_irreps'])
            for imo, sym in enumerate(nos.syms):
                for iirr in self.ioptions['ignore_irreps']:
                    if iirr in sym:
                        #print "%s -> 0"%sym
                        nos.occs[imo] = 0.

        if not(self.ioptions['min_bf'] == ()):
            ntake = 0
            ibf, vbf = self.ioptions['min_bf']
            print("Selecting only NOs with %s > %f ..."%(nos.bf_labels[ibf], vbf))
            for imo in range(nos.ret_num_mo()):
                pop = nos.ret_mo_pop(imo, dosum=2)
                if pop[ibf] > vbf:
                    if self.ioptions['lvprt'] >= 2: print(("Selecting MO %i, occ = %.3f"%(imo+1, nos.occs[imo])))
                    ntake += 1
                else:
                    nos.occs[imo] = 0.
                    if self.ioptions['lvprt'] >= 2: print(("Discarding MO %i"%(imo+1)))
            print("%i NOs selected."%ntake)

        T = numpy.dot(ref_mos.ret_mo_mat(trnsp=False, inv=True), nos.mo_mat)

        state['sden'] = numpy.dot(T,
                                  numpy.dot(numpy.diag(nos.occs), T.transpose()))

        if self.ioptions['unpaired_ana']:
            nu_list = [min(occ, 2.-occ) for occ in nos.occs]
            state['nu'] = sum(nu_list)
            state['nu_den'] = numpy.dot(T,
                                  numpy.dot(numpy.diag(nu_list), T.transpose()))

            nunl_list = [occ*occ*(2-occ)*(2-occ) for occ in nos.occs]
            state['nunl'] = sum(nunl_list)
            state['nunl_den'] = numpy.dot(T,
                                  numpy.dot(numpy.diag(nunl_list), T.transpose()))

class file_parser_rassi(file_parser_libwfa):
    def read(self, mos):

        (energies, oscs, tmp_state_list) = self.read_rassi_output(self.ioptions['rfile'])

    #    if not (self.ioptions['read_libwfa'] or os.path.exists('TRD')):
    #        errmsg= """read_libwfa=False and no TRD directory present!
    #To run the job you have to:
    #1. Run &RASSI specifying the TRD1 keyword
    #2. Run &WFA
    #3. set read_libwfa=True in dens_ana.in
    #
    #or
    #1. Run &RASSI specifying the TRD1 keyword
    #2. call 'mkdir TRD && cp $WorkDir/TRD2* TRD'"""
    #
    #        raise error_handler.MsgError(errmsg)

        if (self.ioptions['read_libwfa']):
            state_list = tmp_state_list
        else:
            print("""Analyzing RASSI job without libwfa.
   WARNING: this only works if all states derive from the same RASSCF calculation!\n""")

            state_list = []
            for lfile in self.ioptions['ana_files']:
                words = lfile.split('_')
                (st1, st2) = (int(words[1]), int(words[2]))

                if self.ioptions['s_or_t'] == 't':
                    if st1 == st2: continue

                    print("Reading %s ..."%lfile)
                    state_list.append({})

                    state_list[-1]['name'] = 'R%i.%i'%(st1, st2)
                    state_list[-1]['exc_en'] = (energies[st1-1] - energies[st2-1]) * units.energy['eV']
                    try:
                        state_list[-1]['osc_str'] = oscs[(st2, st1)]
                    except:
                        print("No osc. strength found for transition %i -> %i"%(st2, st1))
                    state_list[-1]['tden'] = self.init_den(mos)

                    self.read_rassi_den(state_list[-1]['tden'], mos, lfile)

                elif self.ioptions['s_or_t'] == 's':
                    if st1 != st2: continue

                    print("Reading %s ..."%lfile)
                    state_list.append({})

                    state_list[-1]['name'] = 'RASSI_%i'%st1
                    state_list[-1]['exc_en'] = (energies[st1-1] - energies[0]) * units.energy['eV']
                    state_list[-1]['sden'] = self.init_den(mos)

                    self.read_rassi_den(state_list[-1]['sden'], mos, lfile, sden=True)
        return state_list

    def read_rassi_output(self, filen):
        """
        Read the standard output of RASSI.
        """
        state_list = []
        energies = []
        oscs = {}
        libwfa = False

        rfile = open(filen,'r')

        while(True):
            try:
                line = next(rfile)
            except:
                break

            # if 'Total energies (spin-free)' in line:
            #     words = rfile.next().split()
            #     while(len(words) > 0):
            #         energies.append(float(words[-1]))
            #         words = rfile.next().split()

            if 'SPIN-FREE ENERGIES' in line:
                line = next(rfile)
                if 'Shifted by' in line: next(rfile)
                next(rfile)
                next(rfile)
                words = next(rfile).split()
                while(len(words) > 0):
                    energies.append(float(words[1]))
                    words = next(rfile).split()

            elif '++ Dipole transition strengths' in line:
                next(rfile)
                next(rfile)
                next(rfile)
                next(rfile)
                next(rfile)
                words = next(rfile).split()
                while(len(words) > 1):
                    if words[1] == 'Max':
                        break
                    ist = int(words[0])
                    jst = int(words[1])
                    osc = float(words[2])
                    oscs[(ist, jst)] = osc
                    words = next(rfile).split()
                if not self.ioptions['read_libwfa']:
                    break

            elif 'RASSI analysis for reference state' in line or 'RASSI analysis for state' in line:
                state_list.append({})
                state = state_list[-1]
                istate = len(state_list)

                state['name'] = line.split()[-1]

                line = next(rfile)
                line = next(rfile)
                state['exc_en'] = float(line.split()[1])
                try:
                    state['osc_str'] = oscs[(1, istate)]
                except:
                    pass

                libwfa = True

            elif 'RASSI analysis for transiton from state' in line or 'RASSI analysis for transition from state' in line: # typo fixed
                words = line.split()
                typ = words[-1][1:-1]
                (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile("%s_ctnum_atomic.om"%typ)
                if abs(exctmp * units.energy['eV'] - state['exc_en']) > 1.e-3:
                    print(" WARNING: inconsistent energies for %s: %.5f/%.5f - skipping."%(typ, exctmp * units.energy['eV'], state['exc_en']))
                else:
                    state['Om'] = om_at.sum()
                    state['OmAt'] = om_at

            elif libwfa:
                self.parse_line(state, line, rfile)

        rfile.close()

        return energies, oscs, state_list

    def read_rassi_den(self, dens, mos, filen, sden=False):
        """
        Read the output of RASSI generated with TRD1.
        No support for symmetry (yet).
        """
        trd1 = False
        imo = -1
        jmo = -1
        rfile = open(filen,'r')
        TRD_string = 'Active TRD1'
        while True:
            try:
                line = next(rfile)
            except StopIteration:
                break

            if 'Multiplicities' in line:
                words = next(rfile).split()
                mult1 = int(words[0])
                mult2 = int(words[1])
                print('Multiplicities: %i, %i'%(mult1, mult2))

                # Read the spin-traced TDM if the multiplicities are the same
                #   and the spin-difference TDM otherwise
                if mult1 == mult2:
                    TRD_string = 'Active TRD1'
                else:
                    TRD_string = 'Active Spin TRD1'
            elif 'Basis functions' in line:
                nbas = int(next(rfile).split()[0])
                assert(nbas == mos.ret_num_mo())
            elif 'Inactive orbitals' in line:
                ninact = int(next(rfile).split()[0])

                if sden:
                    for imo in range(ninact): dens[imo, imo] = 2.0
                imo = ninact
                jmo = ninact
            elif 'Active orbitals' in line:
                nact = int(next(rfile).split()[0])
            elif TRD_string in line:
                line = next(rfile)
                trd1 = True
            elif trd1:
                strl = len(line) // 5
                for i in range(5):
                    val = float(line[i*strl:(i+1)*strl].replace('D','E'))
                    dens[jmo, imo] = val
                    if 0.2 < abs(val) < 1.9999: print("(i,j)=(%2i,%2i), val=%6.3f"%(imo,jmo,val))

                    jmo += 1
                    if jmo == ninact + nact:
                        jmo = ninact
                        imo += 1

                        if imo == ninact + nact:
                            #print "Finished parsing"
                            trd1 = False
                            break
        rfile.close()
        if trd1:
            raise error_handler.MsgError('Parsing of RASSI output')
        if not sden:
            dens *= 2**(-.5)
