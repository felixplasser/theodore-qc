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

    def print_data(self):
        for NICS_point in self.NICS_data:
            print(NICS_point)

class NICS_point:
    """
    Container for NICS_data.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.NICS_iso = None
        self.NICS_tensor = None

    def __str__(self):
        outstr = "NICS value at (%.5f, %.5f, %.5f)"%(self.x, self.y, self.z)
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
        self.NICS_tensor = NICS_tensor

class NICS_parser_g09(NICS_parser):
    """
    Parse Gaussian NICS calculations.
    """
    def read(self, logfile):
        with open(logfile, 'r') as f:
            Bqind = 0
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    print("Finished parsing %s"%logfile)
                    break

                if line[0:4] == ' Bq ':
                    (x, y, z) = map(float, line.split()[1:4])
                    self.NICS_data.append(NICS_point(x, y, z))
                    
                if 'Bq   Isotropic' in line:
                    NICS_iso = float(line.split()[4])
                    self.NICS_data[Bqind].set_iso(NICS_iso)
                    Bqind += 1