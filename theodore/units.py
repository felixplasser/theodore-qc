"""
version 1.0
author: Felix Plasser
Library containing unit conversion factors. Magnitudes correspond to 1 atomic unit, arranged in dictionaries.
For example, 'eV':27.2113961 can be understood as 27.211 eV/au
"""

from __future__ import print_function, division

energy={'au':1.0,
        'eV':27.2113961,
        'kJ/mol':2625.49963,
        'J':4.359744E-18,
        'kcal/mol':627.50947,
        'rcm':219474.631,
        'nm':45.5633515,
        'Hz':6.579684E15,
        'K':3.157746E+05,
        'RT':3.157746E+05/298.
}
length={'au':1.0,
        'A':0.529177249,
        'm':0.529177249E-10,
        'cm':0.529177249E-8
}
time={'au':1.0,
      'fs':0.0241888432,
      's':0.0241888432E-15
}

mass={'au':1.0,
      'kg':9.10938188E-31,
      'amu':1822.889 # atomic mass unit in a.u. = 1/constants['Nl']/mass['kg']/1000
}
dipole={'D':2.54174619,
        'Cm':8.47835267E-30
}
constants={'Nl':6.02214179E23,
           'c0':137.035999 # speed of light in a.u., 1/alpha
}

# derived quantities
tpa={}
tpa['GM']=10**50 * length['cm']**4 * time['s']

# shortcuts

def eV2nm(val):
    """
    Convert eV to nanometers.
    """
    return energy['eV']*energy['nm']/val

def nm2eV(val):
    """
    Convert nanometers to eV.
    """
    return eV2nm(val)

def eVdiff(au1,au2):
    """
    Energy difference in eV.
    """
    print("%.5f eV"%((au2-au1)*energy['eV']))

def print_units(pformat='%15s : %E'):
    """
    Print information about all the units.
    """
    print(' *** Conversion from atomic units ***')
    print('  Energy')
    for att, val in energy.items():
        print(pformat%(att, val))
    print('  Length')
    for att, val in length.items():
        print(pformat%(att, val))
    print('  Time')
    for att, val in time.items():
        print(pformat%(att, val))
    print('  Mass')
    for att, val in mass.items():
        print(pformat%(att, val))
    print('  Dipole')
    for att, val in dipole.items():
        print(pformat%(att, val))
    print('  Constants')
    for att, val in constants.items():
        print(pformat%(att, val))
    print('  Two-photon absorption')
    for att, val in tpa.items():
        print(pformat%(att, val))

class converter:
    def __init__(self):
        for att, val in energy.items():
            setattr(self, att, val)
        for att, val in length.items():
            setattr(self, att, val)
        for att, val in time.items():
            setattr(self, att, val)
        for att, val in mass.items():
            setattr(self, att, val)
        for att, val in dipole.items():
            setattr(self, att, val)
        for att, val in constants.items():
            setattr(self, att, val)
        for att, val in tpa.items():
            setattr(self, att, val)

u = converter()

if __name__=='__main__':
    print_units()
