import os
from .actions import Action
from .theotools import timeit
from .. import theo_header, lib_tden, lib_exciton, input_options


class AnalyzeTden(Action):
    """
    *** This is the docstring for analyze_tden. ***
    analyze_tden is used for analyzing transition density matrices.
    - Charge-transfer numbers
    - Natural transition orbitals
    - Exciton sizes
    """

    name = 'analyze_tden'

    _colt_description = 'Transition density matrix analysis'

    _questions = """
    ifile = dens_ana.in :: existing_file, alias=f
    """

    @timeit
    def run(ifile):
        ioptions = input_options.tden_ana_options(ifile)
        theo_header.print_header('Transition density matrix analysis', ioptions=ioptions, cfile=__name__)

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

class AnalyzeTdenUnr(Action):

    name = 'analyze_tden_unr'

    _colt_description = 'Transition density matrix analysis'

    _questions = """
    ifile = dens_ana.in :: existing_file, alias=f
    """

    @timeit
    def run(ifile):
        ioptions = input_options.tden_ana_options(ifile)
        theo_header.print_header('Transition density matrix analysis (UHF/UKS)', ioptions=ioptions, cfile=__name__)

        ioptions['jmol_orbitals'] = False

        # ALPHA spin
        print("\nRunning alpha-spin analysis in directory ALPHA")
        ioptions['spin'] = 1

        tdena_alpha = lib_tden.tden_ana(ioptions)
        if 'mo_file' in ioptions: tdena_alpha.read_mos(spin=1)
        tdena_alpha.read_dens()

        try:
            os.mkdir('ALPHA')
        except OSError:
            pass
        os.chdir('ALPHA')
        if 'at_lists' in ioptions:
            tdena_alpha.compute_all_OmFrag()
            if ioptions['print_OmFrag']: tdena_alpha.fprint_OmFrag()
        if ioptions['comp_ntos']:  tdena_alpha.compute_all_NTO()

        print("\n *** ALPHA-spin results ***")
        tdena_alpha.print_summary()
        os.chdir('..')

        # BETA spin
        print("\nRunning beta-spin analysis in directory BETA")
        ioptions['spin'] = -1

        tdena_beta = lib_tden.tden_ana(ioptions)
        if 'mo_file' in ioptions: tdena_beta.read_mos(spin=-1)
        tdena_beta.read_dens()

        try:
            os.mkdir('BETA')
        except OSError:
            pass
        os.chdir('BETA')
        if 'at_lists' in ioptions:
            tdena_beta.compute_all_OmFrag()
            if ioptions['print_OmFrag']: tdena_beta.fprint_OmFrag()
        if ioptions['comp_ntos']:  tdena_beta.compute_all_NTO()

        print("\n *** BETA-spin results ***")
        tdena_beta.print_summary()
        os.chdir('..')

        # ALPHA+BETA
        print("Starting spin-summed analysis")
        # Add the alpha values on top of the beta values
        for i, state in enumerate(tdena_beta.state_list):
            # Add the things that are additive
            for aprop in ['Om', 'OmAt', 'OmFrag', 'S_HE']:
                if aprop in state:
                    state[aprop] += tdena_alpha.state_list[i][aprop]
            try:
                state['Z_HE'] = 2.**(state['S_HE'])
            except KeyError:
                pass

            # Delete the things that are non-additive
            for dprop in ['tden', 'PRNTO', 'Om_desc']:
                if dprop in state:
                    del state[dprop]

        if 'at_lists' in ioptions:
            tdena_beta.compute_all_OmFrag()
            if ioptions['print_OmFrag']: tdena_beta.fprint_OmFrag()

        print("\n *** Spin-summed results ***")
        tdena_beta.print_summary()
