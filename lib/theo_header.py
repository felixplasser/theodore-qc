width=80
    
def print_header(*args):
    print ret_header(*args)
    
def ret_header(title=None):
    hstr  = width*'=' + '\n'
    
    hstr += addlinec("TheoDORE 1.1.4")
    hstr += addlinec("Theoretical Density, Orbital Relaxation and Exciton analysis")
    hstr += addlinec("Felix Plasser")
    
    hstr += width*'-' + '\n'
    
    hstr += addlinec("References")
    hstr += addlinel("Transition density matrix analysis:", 3)
    hstr += addlinel("F. Plasser and H. Lischka")
    hstr += addlinel("J. Chem. Theo. Comp. (2012), 8, 2777.")
    hstr += addlinec()
    hstr += addlinel("Transition and difference density matrix analysis:", 3)
    hstr += addlinel("F. Plasser, S. A. Baeppler, M. Wormit, A. Dreuw")
    hstr += addlinel("J. Chem. Phys. (2014), 141, 024106;")
    hstr += addlinel("J. Chem. Phys. (2014), 141, 024107.")
    hstr += addlinec()
    hstr += addlinel("Exciton analysis:", 3)
    hstr += addlinel("S. A. Baeppler, F. Plasser, M. Wormit, A. Dreuw")
    hstr += addlinel("Phys. Rev. A (2014), 90, 052521.")
    hstr += addlinec()
    hstr += addlinel("Program citation:", 3)
    hstr += addlinel("F. Plasser \"TheoDORE: a package for theoretical density, orbital")
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
