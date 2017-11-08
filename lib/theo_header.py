width=80
    
def print_header(*args, **kwargs):
    print(ret_header(*args, **kwargs))
    
def ret_header(title=None, ioptions=None, ver='1.5.1'):
    hstr  = width*'=' + '\n'
    
    hstr += addlinec("TheoDORE %s"%ver)
    hstr += addlinec("Theoretical Density, Orbital Relaxation and Exciton analysis")
    hstr += addlinec("Felix Plasser")
    
    hstr += width*'-' + '\n'
    
    hstr += addlinec("References for the modules used")
    hstr += addlinec("(see also http://theodore-qc.sourceforge.net/literature.html)")
    hstr += addlinec()
    hstr += addlinel("Transition density matrix analysis:", 3)
    hstr += addlinel("F. Plasser and H. Lischka")
    hstr += addlinel("J. Chem. Theory Comput. (2012), 8, 2777.")
    hstr += addlinec()
    hstr += addlinel("Transition and difference density matrix analysis:", 3)
    hstr += addlinel("F. Plasser, M. Wormit, A. Dreuw")
    hstr += addlinel("J. Chem. Phys. (2014), 141, 024106.")
#    hstr += addlinel("F. Plasser, S. A. Baeppler, M. Wormit, A. Dreuw")
#    hstr += addlinel("J. Chem. Phys. (2014), 141, 024107.")

    hstr += add_exciton(ioptions)
    hstr += add_entanglement(ioptions)
    hstr += add_cclib(ioptions)

    hstr += addlinec()
    hstr += addlinel("Program citation:", 3)
    hstr += addlinel("F. Plasser \"TheoDORE %s: a package for theoretical density, orbital"%ver)
    hstr += addlinel("relaxation, and exciton analysis\"; available from")
    hstr += addlinel("http://theodore-qc.sourceforge.net")

    if not title==None:
        hstr += width*'-' + '\n'
        hstr += addlinec(title)
        
    
    hstr += width*'=' + '\n'
    
    return hstr

def addlinec(line=""):
    return "|" + line.center(width-2) + "|\n"

def addlinel(line="", lpad=5):
    return "|" + lpad*' ' + line.ljust(width-2-lpad) + "|\n"

def add_exciton(ioptions):
    try:
        prop_list = ioptions['prop_list']
    except TypeError:
        return ''

    rstr = ''

    if ('RMSeh' in prop_list) or ('dexc' in prop_list):
        rstr += addlinec()
        rstr += addlinel("Exciton analysis:", 3)
        rstr += addlinel("S. A. Baeppler, F. Plasser, M. Wormit, A. Dreuw")
        rstr += addlinel("Phys. Rev. A (2014), 90, 052521.")

    if 'RMSeh' in prop_list:
        rstr += addlinec()
        rstr += addlinel("Approximate RMSeh/dexc formula:", 3)
        rstr += addlinel("S. A. Mewes, J.-M. Mewes, A. Dreuw, F. Plasser")
        rstr += addlinel("Phys. Chem. Chem. Phys. (2016), 18, 2548.")

    if 'dH-E' in prop_list or 'Corr' in prop_list or 'sigH' in prop_list or 'sigE' in prop_list:
        rstr += addlinec()
        rstr += addlinel("Statistical analysis of excitations:", 3)
        rstr += addlinel("F. Plasser, B. Thomitzni, S. A. Baeppler et al.")
        rstr += addlinel("J. Comp. Chem. (2015), 36, 1609.")

    return rstr

def add_entanglement(ioptions):
    try:
        prop_list = ioptions['prop_list']
    except TypeError:
        return ''

    rstr = ''

    if ('S_HE' in prop_list) or ('Z_HE' in prop_list):
        rstr += addlinec()
        rstr += addlinel("Electron-hole entanglement:", 3)
        rstr += addlinel("F. Plasser")
        rstr += addlinel("J. Chem. Phys. (2016), 144, 194107.")

    return rstr

def add_cclib(ioptions):
    try:
        rtype = ioptions['rtype'].lower()
    except TypeError:
        return ''

    rstr = ''

    if rtype in ['cclib', 'gamess', 'orca']:
        rstr += addlinec()
        rstr += addlinel("cclib for structure parsing:", 3)
        rstr += addlinel("N. M. O'Boyle, A. L. Tenderholt, K. M. Langner")
        rstr += addlinel("J. Comp. Chem. (2008), 29, 839.")

    return rstr
