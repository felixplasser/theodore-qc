from __future__ import print_function, division

from . import lib_struc, error_handler, units
import numpy

class exciton_analysis:
    """
    Perform analysis of an effective exciton wavefunction.
    Approximate atom centered solutions.
    """
    # TODO: add correlation coefficient and covariance
    #   use the same moment construction as in the JCC paper, only discretized
    
    def __init__(self):
        self.distmat = None
    
    def get_distance_matrix(self, struc):
        self.distmat = struc.ret_distance_matrix()
                        
    def ret_RMSeh(self, Om, OmAt):
        """
        Return the root mean square electron-hole distance (Ang).
        """
        if not type(self.distmat) is numpy.ndarray:
            raise error_handler.MsgError("Compute the distance matrix first!")
        
        try:
            MS_dist = numpy.dot(OmAt.flatten(), self.distmat.flatten()**2.) / Om
        except:
            print("\n Error when computing MS_dist!")
            print(" Please check the coordinate file.")
            print(" OmAt: %i x %i"%(len(OmAt), len(OmAt[0])))
            print(" distmat: %i x %i"%(len(self.distmat), len(self.distmat[0])))
            raise
        
        RMS_dist = numpy.sqrt(MS_dist)
        
        return RMS_dist

    def ret_MAeh(self, Om, OmAt):
        """
        Return the mean absolute electron-hole distance (Ang).
        """
        if not type(self.distmat) is numpy.ndarray:
            raise error_handler.MsgError("Compute the distance matrix first!")
        
        MA_dist = numpy.dot(OmAt.flatten(), self.distmat.flatten()) / Om
        
        return MA_dist

    def ret_Eb(self, Om, OmAt, Eb_diag=1.0):
        """
        Return an approximate exciton binding energy (eV).
        """
        if not type(self.distmat) is numpy.ndarray:
            raise error_handler.MsgError("Compute the distance matrix first!")
        
        Eb_dist = self.distmat.flatten() / units.length['A']
        
        for i in range(len(self.distmat)):
            Eb_dist[i + i*len(self.distmat)] = Eb_diag
            
        Eb_au = numpy.dot(OmAt.flatten(), Eb_dist**-1.) / Om
        
        return Eb_au * units.energy['eV']
