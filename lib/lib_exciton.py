import lib_struc
import numpy

class exciton_analysis:
    """
    Perform analysis of an effective exciton wavefunction.
    Approximate atom centered solutions.
    """
    # TODO: add correlation coefficient and covariance (?)
    
    def __init__(self):
        self.distmat = None
    
    def get_distance_matrix(self, coor_file, coor_type):
        struc = lib_struc.structure()
        struc.read_file(coor_file, coor_type)
        self.distmat = struc.ret_distance_matrix()
                
    def ret_RMSeh(self, Om, OmAt):
        """
        Return the root mean square electron-hole distance (Ang).
        """
        if self.distmat == None:
            print " ERROR: compute the distance matrix first!"
            exit(4)
        
        MS_dist = numpy.dot(OmAt.flatten(), self.distmat.flatten()**2.) / Om
        
        RMS_dist = numpy.sqrt(MS_dist)
        
        return RMS_dist
