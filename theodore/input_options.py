"""
Utilities for reading and writing options from/to an input file.
"""

from __future__ import print_function, division
import sys
from . import error_handler


def user_input_py2(inpstr):
    return raw_input(inpstr)


def user_input_py3(inpstr):
    return input(inpstr)


if sys.version_info[0] == 2:
    user_input = user_input_py2
else:
    user_input = user_input_py3


class options:
    """
    Base class for handling input options.
    """
    def __init__(self, ifile):
        self.opt_dict = {}
        self.ifile = ifile

    def __getitem__(self, option):
        return self.get(option, strict=True)

    def get(self, option, strict=True):
        """
        Return the value of an option.
        """
        self.chk_option(option)

        if strict and self.opt_dict[option] == None:
            raise error_handler.MsgError('Option "%s" not defined in file %s!'%(option, self.ifile))
        else:
            return self.opt_dict[option]

    def __setitem__(self, key, val):
        self.opt_dict[key] = val

    def __contains__(self, option):
        """
        Check if an option has been set.

        Raise an error if the option does not even exist.
        """
        self.chk_option(option)

        return self.opt_dict[option] != None

    def has_key(self, option):
        return self.__contains__(option)

    def chk_option(self, option):
        if not option in self.opt_dict:
            raise error_handler.MsgError("Option %s not known!"%option)

    def check_at_lists(self, at_lists, prt_lvl=0):
        """
        Check if an at_lists definition of molecular fragments is useful.
        """
        num_lists = len(at_lists)

        lens = []
        sum_list = []
        for at_list in at_lists:
            sum_list+=at_list
            lens.append(len(at_list))

        numen = len(sum_list)
        maxen = max(sum_list)

        if prt_lvl >= 1:
            print('\nChecking whether the at_lists definition is valid ...')
            if prt_lvl >= 2:
                print('at_lists=', at_lists)
            print('  %i lists with individual numbers of entries:'%(num_lists))
            print(lens)

            print('  %i total entries, with maximal value %i'%(numen,maxen))

        for i in range(1,maxen+1):
            ci = sum_list.count(i)
            if ci!=1:
                print(' WARNING: value %i present %i times in at_lists!'%(i,ci))

    def copy(self, coptions):
        """
        Copy information from a different options instance.
        """
        for key, val in coptions.opt_dict.items():
            self[key] = val

class read_options(options):
    """
    General class for handling input options read from file.
    """
    def __init__(self, ifile, check_init=True):
        options.__init__(self, ifile)

        self.set_defaults()
        self.init = self.read_ifile()

        if check_init: self.check_init()

        self.post_process()

    def check_init(self):
        """
        Check if the instance was properly initialized (the file was read).
        """
        if self.init > 0:
            print("\n ERROR: Input file %s not found!"%self.ifile)
            print("  Please create this file using theoinp")
            exit(0)

    def set_defaults(self):
        """
        Set defaults for the options.
        All possible options should appear here.
        -> inherit for specific implementations
        """
        pass

    def read_ifile(self):
        """
        Read the input file self.ifile.
        Key and value are separated by '='.
        Leading and trailing whitespace is removed.
        """
        try:
            fileh = open(self.ifile, 'r')
        except:
            return 1

        for line in fileh:
            # take out possible comments
            if '#' in line: continue

            words = line.strip().split('=')
            if len(line.strip()) == 0: continue

            if len(words) != 2:
                print(" ERROR: in file %s\n   line cannot be parsed:"%self.ifile)
                print(len(line))
                print(line)
                exit(6)

            key = words[0].strip()

            if words[1] == '':
                raise error_handler.MsgError('Please specify a value for "%s=" in %s!'%(key, self.ifile))

            val = eval(words[1])

            # every possible option has to be initiliazed in set_defaults to avoid confusion
            if not key in self.opt_dict:
                raise error_handler.MsgError('Unknown option in %s: %s'%(self.ifile, key))

            self.opt_dict[key] = val

        return 0

    def get_def(self, option, default):
        self.chk_option(option)

        if self.opt_dict[option] == None:
            return default
        else:
            return self.opt_dict[option]

    def post_process(self):
        pass

class write_options(options):
    """
    General class for writing options to an input file.
    """
    def __init__(self, ifile):
        options.__init__(self, ifile)

        self.ostr = ''

    def read_str(self, title, key, *args, **kwargs):
        """
        Read a string from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_str(titlek, *args, **kwargs)

        self.write_option(key, val)

    def ret_str(self, title, default='', autocomp=False):
        # readline, which is used for auto completion,
        #   creates weird output of the form [?1034h
        #   it should only be imported here
        import readline

        print()
        print(title)

        acstr = ' (autocomplete enabled)' if autocomp else ''
        inpstr = 'Choice%s: '%acstr
        if not default=='': inpstr += '[%s] '%default

        if autocomp:
            readline.set_completer_delims(' \t\n;')
            readline.parse_and_bind("tab: complete")    # activate autocomplete
        val = user_input(inpstr)
        readline.parse_and_bind("tab: ")            # deactivate autocomplete

        if val=='': val = default

        return val

    def read_float(self, title, key, default=1.111):
        """
        Read a float from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_float(titlek, default)

        self.write_option(key, val)

    def ret_float(self, title, default=1.111):
        print()
        print(title)

        inpstr = 'Choice: '
        if not default==1.111: inpstr += '[%f] '%default

        sval = user_input(inpstr)
        if sval=='':
            val = default
        else:
            val = float(sval)

        return val

    def read_int(self, title, key, idef=-1):
        """
        Read a string from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_int(titlek, idef)

        self.write_option(key, val)

    def ret_int(self, title, idef=-1):
        print()
        print(title)

        return self.inp_int(idef)

    def inp_int(self, idef=-1):
        inpstr = 'Choice: '
        if not idef==-1:
            inpstr += '[%i] '%idef

        retval = idef
        while True:
            try:
                retval = int(user_input(inpstr))
            except:
                if retval==-1:
                    print("Please enter an integer number!")
            if retval!=-1: break

        return retval

    def read_yn(self, title, key, default=False):
        """
        Read Boolean from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_yn(titlek, default)

        self.write_option(key, val)

        return val

    def ret_yn(self, question, default=False):
        """
        Ask a yes/no question and return True or False.
        """
        print()
        print(question)

        inpstr = 'Choice (y/n): '
        if default:
            inpstr += '[y] '
        else:
            inpstr += '[n] '

        answer = user_input(inpstr)

        if default:
            return not 'n' in answer.lower()
        else:
            return 'y' in answer.lower()

    def choose_list(self, title, key, opt_expl, default=''):
        """
        Choose an option from a list containing options and explanations.
        """
        titlek = "%s (%s):"%(title, key)

        expl = ["%10s - %s"%(opt, expl) for opt, expl in opt_expl]

        idef = -1
        for ioe, oe in enumerate(opt_expl):
            if oe[0] == default:
                idef = ioe + 1
        ichoice = self.ret_choose_list(titlek, expl, idef)

        val = opt_expl[ichoice-1][0]

        self.write_option(key, val)

    def ret_choose_list(self, title, expl, idef=-1):
        """
        Choose an option from a list containing explanations and return the answer.
        """
        print()
        print(title)

        self.print_list(expl)

        return self.inp_int(idef)

    def print_list(self, plist):
        """
        Print an indexed list to screen.
        """
        iopt = 0
        for p in plist:
            iopt += 1
            print("  [%2i] %s"%(iopt, p))

    def write_list(self, key, wlist, lformat="%i"):
        # write_option can be called directly
        self.write_option(key, wlist)

    def write_option(self, key, val):
        self[key] = val

        if type(val) is str:
            self.ostr += "%s='%s'\n"%(key, str(val))
        else:
            self.ostr += "%s=%s\n"%(key, str(val))

    def flush(self, lvprt=0, choose_file=False):
        if choose_file:
            act_ifile = self.ret_str('Name of input file', self.ifile)
        else:
            act_ifile = self.ifile

        fileh = open(act_ifile, 'w')
        fileh.write(self.ostr)
        fileh.close()
        if lvprt==1:
            print('Finished: File %s written.'%act_ifile)

class dens_ana_options(read_options):
    """
    Input options for density analysis.
    """
    def set_defaults(self):
        # General options
        self['lvprt'] = 1 # print level

        # Read options
        self['mo_file'] = None
        self['rtype']   = None # type of input
        self['rfile']   = None # file to read
        self['ana_files'] = [] # list of files to analyze
        self['read_binary'] = False # read binary files rather than standard output (if applicable)
        self['read_libwfa'] = False # switch to libwfa output (applicable for qctddft, rassi)
        self['s_or_t'] = None # state or transition density matrix analysis
        self['ignore_irreps'] = [] # ignore irreps in the MO file
        self['min_bf'] = () # minimal contribution of a basis function type in the MO file, e.g. (2, 0.5)
        self['rd_ene'] = False # interpret energies as occupations in the NO files
        self['occ_fac'] = 1. # Multiply NO occpuations by this factor
        self['unrestricted'] = False # Read unrestricted orbitals
        self['ana_states'] = [] # Analyze only a set of states (list starts with 1)

        # Output options
        self['output_file']   = "ana_summ.txt"
        self['jmol_orbitals'] = True  # output orbitals in jmol format?
        self['molden_orbitals'] = False  # output orbitals in molden format?
        self['min_occ'] = 0.05 # Minimal occupation for orbital print out
        self['alphabeta'] = False # use alpha/beta rather than neg./pos. to code for hole/electron?
        self['mcfmt']          = '% 10E' # format for molden coefficients
        self['output_prec']   = (7,3) # number of digits and decimal digits for output summary
        self['print_sorted']  = True  # final output sorted by energies

        # tden analysis
        self['Om_formula'] = 1
        self['prop_list'] = []
        self['print_OmAt'] = False   # print the atomic Omega matrix to an .npy file and use for automatic restart
        self['print_OmFrag'] = True # print out the fragment Omega matrix to an ASCII file
        self['eh_pop'] = 1 # print e/h populations: 1 - for fragments, 2 - also for atoms
        self['comp_ntos'] = True
        self['comp_dntos'] = False # Compute the domain NTOs
        self['dnto_frags'] = [] # Compute DNTOs only for these fragments

        # sden analysis
        self['pop_ana'] = True
        self['unpaired_ana'] = True
        self['NO_ana'] = True
        self['AD_ana'] = True
        self['BO_ana'] = True
        self['min_BO'] = 0.5 # minimal bond order to print
        self['mo_pop_type'] = -1

        # options for orbkit
        self['cube_orbitals'] = False  # output orbitals as cube files?
        self['vmd_ntos'] = False
        self['vmd_ntos_iv'] = 0.01
        self['comp_p_h_dens'] = False
        self['vmd_ph_dens'] = False
        self['vmd_ph_dens_iv'] = 0.01
        self['comp_rho0n'] = False
        self['vmd_rho0n'] = False
        self['vmd_rho0n_iv'] = 0.01
        self['comp_rho'] = False # compute densities and unpaired densities
        self['numproc'] = 1
        self['comp_dnto_dens'] = 0 # compute cube files for DNTO densities
            # 0 - none, 1 - only hole, 2 - only electron, 3 - both
        self['fchk_dnto_dens'] = 0  # Print densities to the fchk file
            # 0 - none, 1 - only hole, 2 - only electron, 3 - both
        self['normalize_dnto_dens'] = False # Normalize the DNTO densities
            # for each fragment

        # Additional information
        # irrep labels for output
        self['irrep_labels'] = ['I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8']
        self['ncore'] = {} # dictionary: number of frozen core orbitals per irrep

        # atoms
        self['at_lists'] = None
        self['frag_lists'] = None
        self['coor_file'] = None
        self['coor_format'] = None

        # Internal switches (these should not be set in the input file)
        self['spin'] = 0 # 0 - RHF orbitals, 1 - alpha spin, -1 - beta spin

        # Program specific options
        self['TDA'] = False


class tden_ana_options(dens_ana_options):
    """
    Input options for transition density analysis.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)
        # Read options
        self['s_or_t'] = 't'

        # Output options
        self['output_file']   = "tden_summ.txt"
        self['prop_list'] = ['Om', 'POS', 'PR', 'CT', 'COH', 'CTnt']

        # exciton analysis options
        self['Eb_diag'] = 1.0

class sden_ana_options(dens_ana_options):
    """
    Input options for state density analysis.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)
        # Read options
        self['s_or_t'] = 's'

        # Output options
        self['output_file']   = "sden_summ.txt"
        self['prop_list'] = ['nu', 'nunl', 'y0', 'y1', 'eta', 'p']

class libwfa_parse_options(dens_ana_options):
    """
    Input for parsing libwfa output.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)

        self['output_file']   = "libwfa_summ.txt"
        self['prop_list'] = ['Om_', 'nu', 'nunl', 'PRNTO', 'Z_HE', 'dexc', 'dH-E', 'sigH', 'sigE', 'COV', 'Corr', 'p', 'sigD', 'sigA']
        self['read_libwfa'] = True

class fcd_ana_options(dens_ana_options):
    """
    Input options for fragment charge difference analysis.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)

        self['state_pair'] = None
