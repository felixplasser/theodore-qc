from .actions import Action
from .theotools import timeit
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_sden = importer.lazy_import_as('..lib_sden', 'lib_sden')
    lib_exciton = importer.lazy_import_as('..lib_exciton', 'lib_exciton')
    input_options = importer.lazy_import_as('..input_options', 'input_options')


class AnalyzeSden(Action):

    name = 'analyze_sden'

    _colt_description = 'State density matrix analysis'

    _user_input = """
    ifile = dens_ana.in :: existing_file, alias=f
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_sden': 'lib_sden',
            '..lib_exciton': 'lib_exciton',
            '..input_options': 'input_options'
    })

    @timeit
    def run(ifile):
        # header
        theo_header.print_header(__class__._colt_description, cfile=__name__)
        #
        ioptions = input_options.sden_ana_options(ifile)

        #--------------------------------------------------------------------------#
        # Parsing and computations
        #--------------------------------------------------------------------------#

        sdena = lib_sden.sden_ana(ioptions)
        if 'mo_file' in ioptions: sdena.read_mos()
        sdena.read_dens()

        if ioptions['NO_ana']:  sdena.compute_all_NO()
        if ioptions['AD_ana']:  sdena.compute_all_AD()
        if ioptions['BO_ana']:  sdena.compute_all_BO()

        if ioptions['comp_rho']: sdena.compute_rho()

        #--------------------------------------------------------------------------#
        # Print out
        #--------------------------------------------------------------------------#
        if ioptions['pop_ana']: sdena.print_all_pop_table()
        if ioptions['BO_ana']:  sdena.print_all_BO()
        if ioptions['mo_pop_type'] > 0: sdena.print_mo_pops(ioptions['mo_pop_type'])

        sdena.print_summary()
