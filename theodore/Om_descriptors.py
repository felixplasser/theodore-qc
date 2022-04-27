"""
Computation and storage of Omega descriptors, which are derived from the charge
transfer numbers.
"""

from __future__ import print_function, division

import numpy

class Om_desc_coll:
    """
    Collection of Omega descriptors.
    """
    def __init__(self, Om, OmFrag):
        self.descriptors = {}
        self.OmFrag = OmFrag
        self.OmNorm = OmFrag / Om
        self.numFrag = len(OmFrag)
        
    def ret_val_string(self, desc_list, oformat=' % 4.3f'):
        ret_str = ''
        
        for desc in desc_list:
            ret_str += oformat%self.ret_desc(desc)
            
        return ret_str
        
    def ret_desc(self, desc):
        """
        Return the value of a descriptor. It is either taken from storage or computed.
        
        Return None if the descriptor is not available.
        """
        if desc in self.descriptors:
            return self.descriptors[desc]
        else:
            return self.compute_desc(desc)                
    
    def compute_desc(self, desc):
        """
        Compute the value of descriptor desc.
        """
        #print "Computing %s ..."%desc
        
        self.descriptors[desc] = 0.
        
        if desc == 'POSi':
            self.descriptors[desc] = \
                sum( \
                    (A+1) * sum(self.OmNorm[A,B] for B in range(self.numFrag)) \
                for A in range(self.numFrag))
                
        elif desc == 'POSf':
            self.descriptors[desc] = \
                sum( \
                    (B+1) * sum(self.OmNorm[A,B] for A in range(self.numFrag)) \
                for B in range(self.numFrag))
                
        elif desc == 'POS':
            self.descriptors[desc] = \
                0.5 * (self.ret_desc('POSi') + self.ret_desc('POSf'))
            
        elif desc == 'CT':
            for A in range(self.numFrag):
                for B in range(A+1, self.numFrag):
                    self.descriptors[desc] += self.OmNorm[A, B] + self.OmNorm[B, A]

        elif desc == 'CT2':
            for A in range(self.numFrag):
                for B in range(A+2, self.numFrag):
                    self.descriptors[desc] += self.OmNorm[A, B] + self.OmNorm[B, A]
                    
        elif desc == 'CTnt':
            self.descriptors[desc] = \
                self.ret_desc('POSf') - self.ret_desc('POSi')
                    
        elif desc == 'PRi':
            self.descriptors[desc] = \
                1. / sum( \
                    sum (self.OmNorm[A,B] for B in range(self.numFrag))**2 \
                for A in range(self.numFrag))
            
        elif desc == 'PRf' or desc == 'EEDL':
            self.descriptors[desc] = \
                1. / sum( \
                    sum (self.OmNorm[A,B] for A in range(self.numFrag))**2 \
                for B in range(self.numFrag))
            
        elif desc == 'PR':
            self.descriptors[desc] = \
                0.5 * (self.ret_desc('PRi') + self.ret_desc('PRf'))
            
        elif desc == 'PRh':
            self.descriptors[desc] = \
                2. / (self.ret_desc('PRi')**-1. + self.ret_desc('PRf')**-1.)

        elif desc == 'DEL':
            self.descriptors[desc] = \
                1. / sum( \
                    sum ((self.OmNorm[A,B] + self.OmNorm[B,A])/2. for B in range(self.numFrag))**2 \
                for A in range(self.numFrag))

        elif desc == 'COH':
            self.descriptors[desc] = \
                1. / sum( \
                    sum(self.OmNorm[A,B]**2. for B in range(self.numFrag)) \
                for A in range(self.numFrag)) / self.ret_desc('PR')
        
        elif desc == 'COHh':
            self.descriptors[desc] = \
                1. / sum( \
                    sum(self.OmNorm[A,B]**2. for B in range(self.numFrag)) \
                for A in range(self.numFrag)) / self.ret_desc('PRh')
            
        elif desc in ['MC', 'LC', 'MLCT', 'LMCT', 'LLCT', 'SIEL']:
            self.compute_trans_met()
            
        else:
            return None
            #print "\n ERROR: descriptor %s not implemented!"%desc
            #exit(7)
            
        return self.descriptors[desc]
            
    def compute_trans_met(self):
        """
        Routines specifically for transition metals.
        """
        self.descriptors['MC'] = self.OmNorm[0, 0]
        
        self.descriptors['LC']   = 0.
        self.descriptors['MLCT'] = 0.
        self.descriptors['LMCT'] = 0.
        self.descriptors['LLCT'] = 0.
        for A in range(1, self.numFrag):
            self.descriptors['LC']   += self.OmNorm[A,A]
            self.descriptors['MLCT'] += self.OmNorm[0,A]
            self.descriptors['LMCT'] += self.OmNorm[A,0]
            for B in range(A+1, self.numFrag):
                self.descriptors['LLCT'] += self.OmNorm[A,B] + self.OmNorm[B,A]
                
        #hpop = numpy.sum(OmFrag, 1)
        if self.numFrag >= 3:
            epop = numpy.sum(self.OmFrag, 0)
            self.descriptors['SIEL'] = -epop[1] + 1./(self.numFrag-2.) * numpy.sum(epop[2:])
        else:
            self.descriptors['SIEL'] = None
