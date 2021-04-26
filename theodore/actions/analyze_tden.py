from .actions import Action
from .theotools import timeit
from .. import theo_header, lib_tden, lib_exciton, input_options


class AnalyzeTden(Action):

    name = 'analyze_tden'

    _questions = """
    f = dens_ana.in :: file
    """

    @timeit
    def run(f):
        ioptions = input_options.tden_ana_options(f)
        theo_header.print_header('Transition density matrix analysis', ioptions=ioptions)

        tdena = lib_tden.tden_ana(ioptions)
        if 'mo_file' in ioptions: tdena.read_mos()

        tdena.read_dens()

        if 'at_lists' in ioptions or ioptions['eh_pop'] >= 1:
            tdena.compute_all_OmAt()

        if 'at_lists' in ioptions:
            tdena.compute_all_OmFrag()
            if ioptions['print_OmFrag']: 
                tdena.fprint_OmFrag()

        if ioptions['comp_ntos']:  
            tdena.compute_all_NTO()
        if ioptions['comp_dntos']: 
            tdena.compute_all_DNTO()
        if ioptions['comp_p_h_dens']: 
            tdena.compute_p_h_dens()
        if ioptions['comp_rho0n']: 
            tdena.compute_rho_0_n()

        if 'RMSeh' in ioptions.get('prop_list') or 'MAeh' in ioptions.get('prop_list') or 'Eb' in ioptions.get('prop_list'):
            exca = lib_exciton.exciton_analysis()
            exca.get_distance_matrix(tdena.struc)
            tdena.analyze_excitons(exca)

        if 'Phe' in ioptions['prop_list']:
            tdena.compute_all_Phe()

        #--------------------------------------------------------------------------#
        # Print-out
        #--------------------------------------------------------------------------#

        tdena.print_all_eh_pop()

        tdena.print_summary()
