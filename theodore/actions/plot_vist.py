from .actions import Action
from .theotools import timeit
from .. import theo_header, lib_NICS, error_handler, cclib_interface, input_options


class PlotVist(Action):

    name = 'plot_vist'

    _questions = """
    vist     = :: ilist
    ofile    = VIST.vmd :: file
    scale    = 1.0 :: float
    do_coor  = false :: bool
    plot_all = True :: bool
    logfiles = :: list
    """

    @timeit
    def run(vist, ofile, scale, do_coor, plot_all):
        theo_header.print_header('Read NICS values and prepare VIST plot', cfile=__file__)
        
        with open(ofile, 'w') as fh:
            pass

        ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
        ioptions['rtype'] = 'cclib'
        
        nv = lib_NICS.NICS_parser_g09()
        for ilog, logfile in enumerate(logfiles):
            nv.read(logfile)
            nv.print_data()
            if do_coor: # create coor file to be read by VMD
                ioptions['rfile'] = logfile
                ccparser = cclib_interface.file_parser_cclib(ioptions)
                struc = cclib_interface.structure_cclib()
                struc.read_cclib(ccparser.data)
                coorf = "coor%i.xyz"%ilog
                struc.make_coord_file(file_path=coorf,file_type='Bqxyz')
                open(ofile, 'a').write("mol new %s\n"%coorf)
            nv.vmd_tensors(ofile, vlist, scale, plot_all)
        
        # Instructions for VMD
        if do_coor:
            print("""
        Input file for VMD and coordinate file(s) created. Now run:
           vmd -e %s
            """%ofile)
        else:
            print("""
        Input file for VMD created. Now do the following:
        1. Open VMD and load coordinate file
        2.   File - Load Visualization State - %s
            """%ofile)
