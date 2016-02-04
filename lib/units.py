"""
version 1.0
author: Felix Plasser
Library containing unit conversion factors. Magnitudes correspond to 1 atomic unit, arranged in dictionaries.
"""

energy={'au':1.0,
        'eV':27.2113961,
        'kJ/mol':2.625500E3,
        'J':4.359744E-18,
        'kcal/mol':627.51,
        'rcm':219474.631,
        'nm':45.5633515,
        'Hz':6.579684E15,
        'K':3.157746E+05,
        'RT':3.157746E+05/298.
}
length={'au':1.0,
        'A':0.529177249,
        'm':0.529177249E-10
}
time={'au':1.0,
      'fs':0.0241888432,
      's':0.0241888432E-12
}

mass={'au':1.0,
      'kg':9.10938188E-31
}
dipole={'D':2.54174619,
        'Cm':8.47835267E-30
}
constants={'Nl':6.02214179E23}

# shortcuts
def eV2nm(val):
    return energy['eV']*energy['nm']/val

def nm2eV(val):
    return eV2nm(val)

def eVdiff(au1,au2):
    print "%.5f eV"%((au2-au1)*energy['eV'])
