"""
Analysis routines in the context of NICS (nucleus independent chemical shift)
calculations.

Author: Felix Plasser
Date: Nov. 2020
"""

import numpy

# Input


# Output

class NICS_parser:
    """
    Parse and visualise NICS calculations.
    """
    def __init__(self, logfile=None):
        self.NICS_data = []
        if not logfile is None:
            self.read(logfile)

    def read(self, logfile):
        """
        Read data from quantum chemistry log file.
        """
        raise NotImplementedError

    def print_data(self, outfile='NICS.txt'):
        with open(outfile, 'w') as f:
            f.write(self.NICS_data[0].get_header())
            for point in self.NICS_data:
                f.write(point.get_data())

    def vmd_tensors(self, filen='NICS.vmd'):
        """
        Draw VMD tensors for all the NICS values parsed.
        """
        af = open(filen, 'w')
        self.af = af
        af.write(\
"""axes location Off
display projection Orthographic
display depthcue off
color Display Background white
menu graphics on
mol modstyle 0 0 Licorice 0.100000 30.000000 30.000000
""")

        for point in self.NICS_data:
            point.vmd_tensor(af)

        af.close()
        print("VMD input written to %s"%af.name)

class NICS_point:
    """
    Container for NICS data for each individual ghost atom.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.orig = numpy.array([x, y, z], float)
        self.NICS_iso = None
        self.NICS_tensor = None

    def __str__(self):
        outstr = "NICS value at (% .5f, % .5f, % .5f)"%(self.x, self.y, self.z)
        if not self.NICS_iso is None:
            outstr += ": %.5f"%self.NICS_iso
        return outstr

    def set_iso(self, NICS_iso):
        """
        Set isotropic NICS value.
        """
        self.NICS_iso = NICS_iso

    def set_tensor(self, NICS_tensor):
        """
        Set NICS tensor.
        """
        self.NICS_tensor = numpy.array(NICS_tensor, float)

    def get_data(self):
        """
        Return data as formatted string.
        """
        outstr = "% 11.7f % 11.7f % 11.7f"%(self.x, self.y, self.z)
        if not self.NICS_iso is None:
            outstr += " % 9.5f"%self.NICS_iso
        if not self.NICS_tensor is None:
            for val in self.NICS_tensor.flatten():
                outstr += " % 9.5f"%val

        return outstr + '\n'

    def get_header(self):
        """
        Return header for output.
        """
        outstr = "%9s  %9s  %9s "%('x', 'y', 'z')
        if not self.NICS_iso is None:
            outstr += "      sigma   "
        if not self.NICS_tensor is None:
            outstr += "   s_XX     s_YX      s_ZX       s_XY      s_YY      s_ZY     s_XZ       s_YZ      s_ZZ"
        return outstr + '\n'

    # VMD part - maybe this should be moved into a separate library
    def vmd_tensor(self, af):
        """
        Create a graphical representation of the NICS tensor to be plotted in VMD.
           Should one take the left or right eigenvectors?
             -> The rigorous way would probably be an SVD to represent the non-symmetric matrix
           Maybe show as tensor glyph -> ellipsoid
        """
        self.af = af

        (evals, coor) = numpy.linalg.eig(self.NICS_tensor)

        if not numpy.isreal(evals).all():
            # If there are complex eigenvalues, then these have to be transformed
            #   Form linear combinations of the associated eigenvectors to get a 2-dim irrep
            print("""
*** WARNING: Transforming complex eigenvectors of NICS tensor
      into a real representation""")
            print("NICS Tensor")
            print(self.NICS_tensor)
            print("Eigenvalues")
            print(evals)
            print("Eigenvectors")
            print(coor)

            compl_inds = []
            real_inds  = []
            for i in range(3):
                if numpy.isreal(evals[i]):
                    real_inds.append(i)
                else:
                    compl_inds.append(i)
            assert(len(compl_inds)==2)
            assert(numpy.conj(evals[compl_inds[0]]) == evals[compl_inds[1]])

            coor_new = numpy.zeros([3,3])
            coor_new[:, real_inds[0]]  = numpy.real(coor[:, real_inds[0]])
            coor_new[:, compl_inds[0]] = 2**(.5) * numpy.real(coor[:, compl_inds[0]])
            coor_new[:, compl_inds[1]] = 2**(.5) * numpy.imag(coor[:, compl_inds[0]])
            coor = coor_new

            print("Real transformation matrix")
            print(coor)
            print("Transformed NICS tensor")
            tmp = numpy.dot(numpy.dot(numpy.linalg.inv(coor), self.NICS_tensor), coor)
            print(tmp)
            evals = numpy.diag(tmp)

        for mu in range(3):
            # The eigenvector is in a column
            #print("A v", mu)
            # Check if the eigenvectors are correct -> yes
              # Av = numpy.dot(self.NICS_tensor, coor[:, mu])
              # print(coor[:, mu], Av, Av / evals[mu])
            #print("Abs")
            #print(abs(coor[:, mu]))

            if self.vmd_color(evals[mu], eps=1.):
                fac = 0.1 * abs(evals[mu])
                self.plot_quad_comp(fac, self.orig, coor[:, mu], 0.1 * fac**.5)
                #self.plot_quad_cone(fac, self.orig, coor[:, mu], 0.5 * fac)

    def plot_quad_comp(self, fac, orig, vec, rad):
        """
        Plot component of quadrupole moment.
        """
        self.af.write('draw cylinder ')
        self.vmd_coors(orig - fac*vec)
        self.vmd_coors(orig + fac*vec)
        self.af.write('radius % .3f\n'%rad)

        self.af.write('draw sphere ')
        self.vmd_coors(orig + fac * vec)
        self.af.write('radius % .3f\n'%(2*rad))

        self.af.write('draw sphere ')
        self.vmd_coors(orig - fac * vec)
        self.af.write('radius % .3f\n'%(2*rad))

    def plot_quad_cone(self, fac, orig, vec, rad):
        """
        Plot component of quadrupole moment as cone.
        """
        self.af.write('draw cone ')
        self.vmd_coors(orig - fac*vec)
        self.vmd_coors(orig)
        self.af.write('radius % .3f\n'%rad)

        self.af.write('draw cone ')
        self.vmd_coors(orig + fac*vec)
        self.vmd_coors(orig)
        self.af.write('radius % .3f\n'%rad)

    def vmd_coors(self, coor):
        self.af.write('{% .3f % .3f % .3f} '%(coor[0], coor[1], coor[2]))

    def vmd_color(self, val, eps=1.e-3):
        if abs(val) < eps:
            return False
        elif val > 0.0:
            self.af.write('draw color blue\n')
            return True
        else:
            self.af.write('draw color red\n')
            return True

class NICS_parser_g09(NICS_parser):
    """
    Parse Gaussian NICS calculations.
    """
    def read(self, logfile, lvprt=1):
        with open(logfile, 'r') as f:
            Bqind = 0
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    print("Finished parsing %s"%logfile)
                    if lvprt >= 1:
                        for point in self.NICS_data:
                            print(point)
                    break

                if line[0:4] == ' Bq ':
                    (x, y, z) = map(float, line.split()[1:4])
                    self.NICS_data.append(NICS_point(x, y, z))
                    
                if 'Bq   Isotropic' in line:
                    NICS_iso = float(line.split()[4])
                    self.NICS_data[Bqind].set_iso(NICS_iso)

                    tensor = []
                    for i in range(3):
                        words = next(f).split()
                        tensor.append([float(words[1]), float(words[3]), float(words[5])])
                    self.NICS_data[Bqind].set_tensor(tensor)

                    Bqind += 1