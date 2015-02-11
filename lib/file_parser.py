"""
Parsing of files produced by different quantum chemical programs.
"""

import units, lib_mo, error_handler
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
    
    def parse_key(self, state, key, line, search_string, ind=-1):
        """
        Find search_string in the specified line and set it as state[key].
        """
        if search_string in line:
            sr_line = line.strip(search_string).replace(',','')
            state[key] = float(sr_line.split()[ind])
            
    def delete_chars(self, line, delete):
        """
        Delete the characters in delete from line.
        """
        retline = line
        for char in delete:
            retline = retline.replace(char, '')
            
        return retline

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
            print "\n ERROR when parsing line:"
            print line
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
               mos.symsort(self.ioptions['irrep_labels'])
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

            state['tden'][iocc,ivirt] = conf.coeff
                
    def set_tden_bin(self, state, mos, lvprt=1):
        """
        Set the transition density matrix elements acoording to the binary files with state information.
        """
        iirrep=self.ioptions['irrep_labels'].index(state['irrep']) + 1
        smult=state['mult']
        istate=state['state_ind']
        #CCfilen = 'CCRE0-%i--%s---%i'%(iirrep, smult, istate)
        CCfilen = 'CCRE0{0:->2d}{1:->3s}{2:->4d}'.format(iirrep, smult, istate)

        if lvprt >= 1: print 'Reading binary file %s ...'%CCfilen

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
        nact = nentry / nvirt
        nfrzc = nocc - nact
        
        if lvprt >= 1:
            print '  Method: %s, number of entries: %i'%(method, nentry)
            print '  num_mo: %i, nocc: %i, nvirt: %i, nact: %i, nfrzc: %i'%(num_mo, nocc, nvirt, nact, nfrzc)
        if lvprt >= 2:
            print '  Syms:', mos.syms

        if nfrzc > 0:
            print '\n\n  WARNING: Frozen core orbitals should be kept out of the molden file when Reading CCRE0* files!\n'
            
        # write the collected data into the correct block of the 1TDM
        for iocc in xrange(nfrzc, nocc):
            for ivirt in xrange(nocc, num_mo):
                state['tden'][iocc, ivirt] = struct.unpack('d', CCfile.read(8))[0]

        lbytes = struct.unpack('4s', CCfile.read(4))
        if lvprt >= 2:
            print '   Last four bytes:', lbytes
        if not CCfile.read(1) == '':
            raise error_hanlder.MSGError('parsing file %s'%CCfilen)

        if lvprt >= 3:
            print 'parsed tden:'
            print state['tden']

    def ret_conf_ricc2(self, rfile='ricc2.out'):
        """
        Return information about configurations in a Turbomole calculation.
        """
        ret_list = []
        curr_osc = 0 # running index for oscillator strengths
        lines = open(rfile,'r') 
        while True: # loop over all lines
          try:
            line = lines.next()
          except StopIteration:
            print "Finished parsing file %s"%rfile  
            break
          else:
            if ' sym | multi | state' in line:
                # parse the multiplicity from here
                multis = []
                line = lines.next()
                line = lines.next()
                line = lines.next()
                while True:
                    line = lines.next()
                    if '---' in line: continue
                    if '===' in line: break
                    
                    words = self.delete_chars(line, ['|']).split()
                    
                    multis.append(words[1])                   
                    
            elif 'Energy:   ' in line:
                ret_list.append({})
                words = line.split()
                ret_list[-1]['exc_en'] = float(words[3])
                
                for i in xrange(3): line = lines.next()
                    
                words = line.split() # type: RE0
                ret_list[-1]['irrep'] = words[4]
                ret_list[-1]['state_ind'] = int(words[6])
                ret_list[-1]['mult'] = multis.pop(0)
                
                for i in xrange(3): line = lines.next()
                    
                ret_list[-1]['char'] = []
                while True: # loop over information about this excitation
                    line = lines.next()
                    
                    if '===' in line: break
                    
                    conf = configuration()
                    conf.read_ricc2_line(line)                    
                    ret_list[-1]['char'].append(conf)
                    
            elif 'oscillator strength (length gauge)' in line:
                words = line.split()
                ret_list[curr_osc]['osc_str'] = float(words[5])
                curr_osc+=1
            elif 'irred. repres.:' in line:
                self.ioptions['irrep_labels'] = line.split()[2:]
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
                        print " Error: invalid occupation!", occ
                        exit(5)

            nocc=len(occmap)
            nvirt=len(virtmap)

            print "\n Considering state: %s"%state['name']
            print "  Number of occupied orbitals:", nocc
            print "  Number of virtual orbitals:", nvirt
            print "  Tensor space dimension:", nocc*nvirt

            if os.path.exists('sing_%s'%state['irrep']):
              readf = 'sing_%s'%state['irrep']
              print '  Reading information of singlet calculation from file %s'%readf
            elif os.path.exists('trip_%s'%state['irrep']):
              readf = 'trip_%s'%state['irrep']
              print '  Reading information of triplet calculation from file %s'%readf
            elif os.path.exists('ciss_%s'%state['irrep']):
              readf = 'ciss_%s'%state['irrep']
              print '  Reading information of TDA singlet calculation from file %s'%readf
            elif os.path.exists('cist_%s'%state['irrep']):
              readf = 'cist_%s'%state['irrep']
              print '  Reading information of TDA triplet calculation from file %s'%readf
            else:
              print 'No file with information about the excited state (sing_a, trip_a, ...) found!'
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
                    words = [line[0+20*i:20+20*i] for i in xrange(4)]                    
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
    def read(self):
        """
        Read the .om files as created with libwfa.
        These already contain OmAt
        """
        state_list = []
        
        basedir='.'
        #ist = 1
        for omfile in os.listdir(basedir):
            suff = omfile.split('.')[-1]
            if suff!='om': continue
           
            (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile(os.path.join(basedir,omfile))
            if typ == None:
                continue

            state_list.append({})
          
            try:
                multl, iirrep, ist = typ.split('_')
            except:
                state_list[-1]['name'] = 'st'
            else:
              
                state_list[-1]['mult'] = {'singlet':'(1)',
                                            'triplet':'(3)',
                                            'any':'(-)'}[multl]                                  
                state_list[-1]['state_ind'] = int(ist)                
                state_list[-1]['irrep'] = self.ioptions.get('irrep_labels')[int(iirrep)]
                
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
        print "Reading: %s ..."%fname
        
        try:
            rfile=open(fname,'r')
        except IOError:
            print "\n WARNING: could not open %s."%fname
            print "Did the calculation converge?\n"
            return None, None, None, None, None, None

        line=rfile.readline()
        words=line.split()
        typ = words[0]
        excen = float(words[1])
        osc = float(words[2])

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
                print " All values read in"
            
        return typ, excen, osc, dima, dimb, outarr
    
class file_parser_qcadc(file_parser_libwfa):
    """
    Parse information from qchem.out in addition to the .om file.
    """
    def read(self):
        state_list = []
        
        basedir='.'
        
        exc_diff = False
        exc_1TDM = False
        for line in open(self.ioptions.get('rfile')):
            words = line.split()
            if 'Irreducible representations in point group:' in line:
                self.irrep_labels = line.lstrip('Irreducible representations in point group:').split()
                
            elif ' Term symbol' in line:
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
                
            elif ' Excitation energy:' in line:
                exc_chk = float(words[2])
                
                if not 'exc_en' in state_list[-1]:
                    state_list[-1]['exc_en'] = exc_chk
                
                if abs(exc_chk - state_list[-1]['exc_en']) > 1.e-4:
                    print exc_chk, state_list[-1]['exc_en']
                    raise error_handler.MsgError("Excitation energies do not match")
                    
            elif 'Exciton analysis of the difference density matrix' in line:
                if len(state_list) > 0: exc_diff = True
                exc_1TDM = False
                
            elif 'Exciton analysis of the transition density matrix' in line:
                exc_diff = False
                if len(state_list) > 0: exc_1TDM = True
                
            elif 'Transition Summary' in line:
                break

            if len(state_list) > 0:
                self.parse_key(state_list[-1], 'dip', line, 'Total dipole')
                self.parse_key(state_list[-1], 'r2', line, 'Total <r^2>')
                self.parse_key(state_list[-1], 'nu', line, 'Number of unpaired electrons:', 2)
                self.parse_key(state_list[-1], 'nunl', line, 'Number of unpaired electrons')                
                self.parse_key(state_list[-1], 'p', line, 'Number of detached / attached electrons')
                self.parse_key(state_list[-1], 'PRNTO', line, 'PR_NTO')
                
            if exc_diff:
                self.parse_key(state_list[-1], 'sigD', line, 'Hole size')
                self.parse_key(state_list[-1], 'sigA', line, 'Electron size')
                self.parse_key(state_list[-1], 'dD-A', line, '|<r_e - r_h>|')

            if exc_1TDM:
                self.parse_key(state_list[-1], 'dexc', line, 'RMS electron-hole separation')
                self.parse_key(state_list[-1], 'sigH', line, 'Hole size')
                self.parse_key(state_list[-1], 'sigE', line, 'Electron size')
                self.parse_key(state_list[-1], 'dH-E', line, '|<r_e - r_h>|')
                self.parse_key(state_list[-1], 'COV', line, 'Covariance(r_h, r_e) [Ang^2]')
                self.parse_key(state_list[-1], 'Corr', line, 'Correlation coefficient')
                
        return state_list

    def om_file_name(self, state):
        """
        Construct the name of the .om file
        """
        multo = {'(1)':'singlet',
                 '(3)':'triplet',
                 '(-)':'any'}[state['mult']]
                 
        irrepo = self.irrep_labels.index(state['irrep'])
        
        if irrepo == 0 and state['mult'] == '(1)':
            # add 1 for ground state irrep
            state_indo = state['state_ind'] - 1
        else:
            state_indo = state['state_ind']
        
        state['fname'] = '%s_%i_%i'%(multo, irrepo, state_indo)
        
        return '%s_%i_%i_ctnum_atomic.om'%(multo, irrepo, state_indo)
    
class file_parser_qctddft(file_parser_base):
    def read(self, mos):
        """
        Read the X vector from standard output. Y is discarded.
        """
        state_list = []
        
        if self.ioptions.get('TDA'):
            ststr = 'TDDFT/TDA Excitation Energies'
        else:
            ststr = 'TDDFT Excitation Energies'
            
        print "Parsing %s for %s ..."%(self.ioptions.get('rfile'), ststr)
        
        tdread = False
        for line in open(self.ioptions.get('rfile'), 'r'):
            if ststr in line:
                tdread = True
            elif 'TDDFT calculation will be performed' in line:
                tdread = False
            elif 'Timing summary' in line:
                tdread = False
                    
            if tdread:                
                words = line.replace(':','').split()
                if 'Excited state' in line:
                    state_list.append({})
                    
                    state_list[-1]['state_num'] = int(words[2])
                    state_list[-1]['exc_en'] = float(words[-1])
                    state_list[-1]['name'] = 'X%i'%state_list[-1]['state_num']
                    
                    state_list[-1]['tden'] = self.init_den(mos, rect=True)
                    
                elif 'Strength' in line:
                    state_list[-1]['osc_str'] = float(words[-1])
                    
                elif 'amplitude' in line:
                    # ignore the Y vector. Which one would be the correct sign?
                    if 'Y:' in line: continue
                    
                    awords = self.delete_chars(line, ['X:', 'Y:', 'D', 'V', '(', ')', '-->']).split()
                    
                    iocc = int(awords[0]) - 1
                    ivirt = int(awords[1]) + mos.ret_ihomo()
                    
                    coeff =  float(awords[4])
                    
                    state_list[-1]['tden'][iocc, ivirt] += coeff
        
        return state_list

#---

class file_parser_col(file_parser_base):
    def read_iwfmt(self,  fname):
        raise error_handler.NIError()
    
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
        if not self.ioptions['ncore'] == 0:
            print """
  WARNING: adjusting MO file for frozen core orbitals!
  This only works without symmetry.
            """
        # change the format of the MO labels
        mos.syms2 = self.ioptions['ncore']*['xxx']+[''.join(sym.lower().split('_')) for sym in mos.syms]
        
        state['tden'] = self.init_den(mos)
        
        tmp = filen.replace('state','').replace('drt','').split('TO')[-1]
        ir_st = tmp.split('.')
        state['irrep'] = self.ioptions.get('irrep_labels')[int(ir_st[0]) - 1]
        state['state_ind'] = int(ir_st[1])
        state['name'] = '%s.%i'%(state['irrep'], state['state_ind'])
        
        rfile = open(filen,'r')
        while True:
            try:
                line=rfile.next()
            except StopIteration:
                break
            if ' symm-trans density matrix:' in line:
                print 'file %s, reading:'%filen
                print line.strip('\n')
                line=rfile.next()
                self.read_block_mat(state, mos, rfile, sym=1)
            elif 'asymm-trans density matrix:' in line:
                print 'reading ...'
                print line[:-1]
                line=rfile.next()
                self.read_block_mat(state, mos, rfile, sym=-1)
            elif 'Transition energy:  ' in line:
                words=rfile.next().split() # use the next line with "eV"
                state['exc_en']=float(words[2])
            elif 'Oscillator strength : ' in line:
                words = line.split()
                state['osc_str']=float(words[-1])
                
    def read_block_mat(self, state, mos, rfile, sym):
        """
        Parse the block matrix output in the listing file.
        This is not the cleanest routine but it seems to work ...
        """
        isqr2 = 1. / numpy.sqrt(2.)
        
        while(1):
            words=rfile.next().replace('MO','').split()
            
            if 'density' in words: # skip two lines in case the next symmetry block is coming
                line=rfile.next()
                words=rfile.next().replace('MO','').split()
            
            try:
                head_inds = [mos.syms2.index(sym2) for sym2 in words]
            except:
                print words
                raise
            
            words=rfile.next().replace('MO','').split()
            while(len(words)>0): # loop over a block with constant header labels
                if 'integral' in words: break
                try:
                    left_ind = mos.syms2.index(words[0])
                except:
                    print '\n ERROR parsing: '
                    print words
                    raise
                    
                for i, word in enumerate(words[1:]):
                    val = float(word)                                                            
                    if abs(val) > 0.1: print "(%2i->%2i)-(%s->%s), val=% 8.4f"%\
                       (head_inds[i],left_ind,mos.syms2[head_inds[i]],mos.syms2[left_ind],val)
                    
                    state['tden'][head_inds[i],left_ind] += isqr2 * val
                    if left_ind != head_inds[i]:
                        state['tden'][left_ind,head_inds[i]] += isqr2 * sym*val
                        
                words=rfile.next().replace('MO','').split()
                
            if 'integral' in words: break
    
class file_parser_col_mrci(file_parser_col):
    def read(self, mos):
        # TODO: separate between sden and tden
        state_list = []
        
        for lfile in os.listdir('LISTINGS'):
            if not 'trncils' in lfile: continue
            
            print "Reading %s ..."%lfile
            state_list.append({})
            self.read_trncils(state_list[-1], mos, 'LISTINGS/%s'%lfile)

        return state_list
    
class file_parser_nos(file_parser_base):
    """
    Interpret the MO-file as a diagonal density.
    """
    def read(self, mos):
        state_list = []
        
        #state_list.append({})
        #self.read_ref_nos(state_list[-1], mos)
        
        for no_file in self.ioptions['no_files']:
            state_list.append({})
            self.read_no_file(state_list[-1], mos, no_file)
            
        for istate, state in enumerate(state_list):
            state['exc_en'] = float(istate + 1) # set fake excitation energy
            state['state_num'] = istate + 1
            
            state['fname'] = self.ioptions['no_files'][istate]
            
            # perform some simplifications to the filename
            # for Columbus:
            tmp_name = state['fname'].replace('MOLDEN/','').replace('molden_no','').rstrip('.sp')
            tmp_name = tmp_name.replace('state', 'S').replace('drt', 'D')
            # for Q-Chem
            tmp_name = tmp_name.replace('NOs/','').replace('.mo', '')
            
            state['name'] = tmp_name

        return state_list    
    
    def read_ref_nos(self, state, ref_nos):
        """
        Read the reference MOs.
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
