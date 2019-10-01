"""
Library with some utilities that do not fit anywhere else.
"""

class cube_file:
    """
    Analyse a cube file and compute isovalues corresponding to volume integrals.
    """
    def __init__(self, fname):
        self.fname = fname
        self.vals =  None
        self.avals = None

    def read(self, lvprt=0):
        if lvprt >= 1:
            print('Analysing %s ...'%self.fname)
        with open(self.fname, 'r') as f:
            line = next(f)
            line = next(f)
            line = next(f)
            natom = int(line.split()[0])
            words = next(f).split()
            x = float(words[1])
            words = next(f).split()
            y = float(words[2])
            words = next(f).split()
            z = float(words[3])
            self.V = x * y * z
            for iat in range(natom):
                line = next(f)

            self.vals = []
            while True: # loop over all lines
                try:
                    line = next(f)
                except StopIteration:
                    if lvprt >= 1:
                        print("Reached end of file %s"%f.name)
                    break

                self.vals += map(float, line.split())

        self.s = sum(self.vals)
        self.avals = [abs(val) for val in self.vals]
        self.abss = sum(self.avals)
        self.sqs = sum(val*val for val in self.vals)
        self.minval = min(self.vals)
        self.maxval = max(self.vals)
        if lvprt >= 1:
            print('Integral: % .6f, Abs. Int.: % .6f, Squ. Int. % .6f:'%(self.s * self.V, self.abss * self.V, self.sqs * self.V))
            print('Min: % .6f, Max % .6f'%(self.minval, self.maxval))
        self.avals.sort(reverse=True)

    def ret_isovals(self, frac=[0.1, 0.5, 0.75, 0.9, 0.95, 0.99], lvprt=0):
        """
        Return the isovalues that correspond to a specific fraction of the density.
        """
        if self.vals is None:
            self.read(lvprt)

        retvals = []

        frac.sort()
        ifrac = 0
        ps = 0.
        thres = frac[ifrac] * self.abss
        for val in self.avals:
            ps += val
            if ps > thres:
                retvals.append(val)
                if lvprt >= 1:
                    print("isoval: %.6f at %.4f (%.6f)"%(val, frac[ifrac], ps/self.abss))
                ifrac += 1
                try:
                    thres = frac[ifrac] * self.abss
                except IndexError:
                    break

        return retvals

    def dot(self, other, lvprt=0):
        """
        Compute the dot product / overlap between two cube files.
        """
        if self.vals is None:
            self.read(lvprt=lvprt)
        if other.vals is None:
            other.read(lvprt=lvprt)
        assert self.V == other.V

        dot = sum(val * other.vals[ival] for ival, val in enumerate(self.vals)) * self.V
        if lvprt >= 1:
            print("Computing dot product between %s and %s"%(self.fname, other.fname))
            print("Dot: % .6f"%dot)

        return dot
