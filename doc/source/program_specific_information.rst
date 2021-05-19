Program specific information 
----------------------------

In the case of TDDFT, "Exact" is used for cases where the whole **X** or **X** + **Y** vector is parsed.

Use the <code>rtype</code> keyword in <code>dens_ana.in</code> to specify the respective method.


In general, TheoDORE needs two types of information: density matrices and MO coefficients. This page contains specific information about how to obtain this information from the different quantum chemistry programs.


Q-Chem
~~~~~~

ADC
___

The ADCMAN module of Q-Chem is already interfaced to the wavefunction analysis library libwfa. Therefore, most analysis steps can be performed within Q-Chem. The main purpose of TheoDORE is to enable [plotting](PLotting) of *electron-hole* correlation plots. Futhermore, specific analysis descriptors can be parsed and put into a convenient table.

To run the libwfa analysis set in the input file:

::

    state_analysis    true
    adc_print    3

TDDFT
_____

The CIS/TDDFT module of Q-Chem is directly interfaced to the wavefunction analysis library libwfa (see [DOI: 10.1063/1.4935178](http://dx.doi.org/10.1063/1.4935178) for more details). Most analysis steps are performed within Q-Chem and the main purpose  of TheoDORE is to enable [plotting](PLotting) of *electron-hole* correlation plots.

To run the libwfa analysis set in the input file:

::

    state_analysis    true

Alternatively, the standard output can be parsed. To do this set <code>read_libwfa=False</code> in <code>dens_ana.in</code>. In this case, the X vector is parsed and interpreted as the 1TDM while the Y vector is ignored. The following non-standard options have to be set in the input file:

::

    cis_ampl_print	1
    print_orbitals	5000
    molden_format	true

The first one assures that enough amplitudes are printed for a meaningful semi-quantitative analysis. The second and third cause the print-out of all orbitals in Molden format, written to the end of the standard output file. This Molden file has to be extracted and copied to a new file <code>qchem.mld</code>.

Per default RPA vectors are parsed. If you are interested in TDA vectors, specify

::

    TDA=True

in <code>dens_ana.in</code>.

fchk files
__________
It is also possible to parse formatted checkpoint (fchk) files generated using the

::

    state_analysis    true

option. In this case, TheoDORE can read the transition density matrices and do the full analysis.

libwfa
______
It is also possible to parse generic <code>libwfa</code> output by using the libwfa keyword.

Columbus
~~~~~~~~

MCSCF
_____
MCSCF transition density matrices are written by the program <code>mcscf.x</code>. To compute these, it is either possible to compute non-adiabatic coupling vectors or to use the shortcut of creating an <code>mcdenin</code> file. This file has the same structure as the <code>transmomin</code> file, e.g.

::

    MCSCF
    1  2  1  1
    1  3  1  1
    1  4  1  1
    1  5  1  1
    1  6  1  1

will compute all transition moments between the 1st and the 2nd through 6th states in DRT 1.

After the computation call <code>write_den.bash</code> to convert the binary files into a form that can be read by TheoDORE. For this purpose the $COLUMBUS variable has to be set.

State density matrix analysis is possible when the computation of dipole moments was requested. In this case it is possible to use the above procedure or to simply analyze the NO coefficient files in the <code>MOLDEN</code> directory.

MR-CI
_____

The MR-CI transition density analysis works in the case that transition moments were requested in the job control stage. After this it is assumed that <code>LISTINGS/trncils....</code> files with the transition density matrices are available. Frozen core orbitals have to be specified explicitly in the input file to get correct results. To specify, for example, four frozen orbitals of *a* symmetry and three of *b* use:

::

    ncore={'a':4, 'b':3}

In addition, the MO-coefficients of the preceeding MCSCF calculation have to be made available, typically in <code>MOLDEN/molden_mo_mc.sp</code>.

For a state density analysis at the MR-CI level, the NO files should be read in rather than using the <code>colmrci</code> functionality. Alternatively, an attachment/detachment analysis can be done with the [densav.x](https://www.univie.ac.at/columbus/docs_COL70/utilities.html#densav) functionality of COLUMBUS.

Molcas
~~~~~~
The prefered way to use Molcas is through the [libwfa library](https://github.com/libwfa/libwfa), which is available through [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas). In this way it is possible to analyze RASSCF and MS-CASPT2 computatoins.

First, run Molcas using the &WFA module

::

    &RASSI
    TRD1

    &WFA
    h5file = $Project.rassi.h5

and copy back the <code>*.om</code> files. In <code>theoinp</code> specify "y" for

::

    Did you use &WFA? (read_libwfa):
    Choice (y/n): [y] 

and proceed as usual.

Molcas (old)
~~~~~~~~~~~~

It is also possible to parse Molcas RASSI output but this only works for singlet states.

RASSCF
______
In the case of Molcas, the output of the RASSI program is parsed. This gives access to RASSCF density and transition density matrices, so far without explicit point group symmetry. **Note**: The output is only parsed correctly if all states in the RASSI computation derive from the same RASSCF computation and if the specified Molden file derives from this calculation.

First run a RASSCF + RASSI job with the (undocumented) <code>TRD1</code> keyword:

::

    &RASSI
        TRD1

Then copy the transition densities to a directory <code>TRD</code>:

::

    mkdir TRD && cp $WorkDir/TRD2* TRD

Alternatively, a state density matrix analysis can be performed by using the natural orbitals created by Molcas. However, for an analysis of unpaired electrons the NOs have to be changed from spin-orbitals to spatial orbitals.

MS-CASPT2
_________

For an MS-CASPT2 calculation, the following input sections can be used

::

    &CASPT2
    multistate = 4 1 2 3 4
    imag = 0.3

    >> SAVE $Project.JobMix JOB001

    &RASSI
    NROFJOBIPHS
    1 4
    1 2 3 4
    CIPR
    TRD1


This will yield density matrices mixed according to the MS-CASPT2 calculation, which can in turn be analyzed by TheoDORE.

Unfortunately, it is not possible to use the EJOB keyword in connection with this procedure. Therefore, the energies and oscillator strengths given are not consistent!

Tubomole
~~~~~~~~

CC2 / ADC(2)
____________

If you have the binary <code>CCRE0*</code> files, written by Turbomole, available, then choose the option

::

    read_binary=True


in <code>dens_ana.in</code>. Use <code>tm2molden</code> without further options to create the MO file. For printing the NTOs, it is not possible to use <code> jmol_orbitals</code> in this case, but only <code>molden_orbitals</code>.

Alternatively, approximate transition density matrices can be read directly from the standard output of <code>ricc2</code>. The MO file is again created with <code>tm2molden</code>. However, it is important that also the frozen orbitals are contained in the MO file. This can be achieved by running the following commands:

::

    #!/bin/bash
    sed -i "/implicit core/d" control
    echo -e "\n\n"|tm2molden

TDDFT
_____
In the TDDFT case, the <code>sing_a</code> or <code>trip_a</code> files are parsed and interpreted as 1TDMs. Unfortunately, this analysis only works if no explicit symmetry is chosen in the initial job setup.

MO-coefficients have to be supplied by <code>tm2molden</code>.

Terachem - TDDFT
~~~~~~~~~~~~~~~~
For a trans. dens. mat. analysis, the CI vectors are read from standard output and the MO coefficients from a Molden file produced by Terachem. To print more CI vector elements, use

::

    cisprintthresh 0.01

A state/difference density matrix analysis is possible by using the NO files produced when using

::

    cisnos   yes

Parsing of NO files
~~~~~~~~~~~~~~~~~~~
NO files can be parsed directly using either <code>analyze_sden.py</code> or <code>analyze_NOs.py</code>. In this case it is important that one reference file is given, which contains the full, invertible MO-matrix.

TheoDORE assumes that the NO files are given with respect to spatial orbitals (occupation between 0 and 2). If spin NOs are given, then the analysis of unpaired electrons will not give suitable results.

ORCA - TDDFT
~~~~~~~~~~~~
Starting in TheoDORE 2.0.1, the preferred version of parsing ORCA TDDFT jobs uses a Molden format file and the  `orca.cis` file.

1. Run an ORCA job and copy back the `orca.gbw` and `orca.cis` files
*Note*: the filename `orca.cis` is hardcoded in TheoDORE
2. Create a molden file using `orca_2mkl orca -molden`
3. Run TheoDORE and select `13` at

::

    Type of job (rtype):
    ...
      [12]      cclib - Use external cclib library: Gaussian, GAMESS, ...
      [13]       orca - ORCA TDDFT (using a Molden file and cclib)
    ...
    Choice: 13

ORCA using cclib
~~~~~~~~~~~~~~~~
Alternatively, ORCA can be parsed entirely with the [cclib library](http://cclib.github.io/). If you want to do that, set the following output options:

::

    %output
     PrintLevel Normal
     Print[ P_MOs ] 1
     Print[ P_Overlap ] 1
    end

It is recommended also in this case to read the CI-vectors from the binary file <code>orca.cis</code> rather than from standard output. To do this, set

::

    read_binary=True

In the case of TDA both options work, for RPA <code>read_binary=True</code> has to be used.

Gaussian - TDDFT
~~~~~~~~~~~~~~~~
Gaussian is parsed with the [cclib library](http://cclib.github.io/). Set the `pop=full iop(9/40=3)` option to increase the number of CI vector elements printed. Use `GFINPUT` to print the basis functions and `iop(3/33=4)` to get the overlap matrix.

Example input:

::

    #p PBEPBE/6-31G* td=(singlets, nstates=10) pop=full iop(9/40=3) GFINPUT

For some applications, in particular in connection with ORBKIT, it is advisable to supply an externally generated molden file with orbital information. For this purpose, open the Gaussian-log file in Molden. Choose "Write - Molden Format" and save as `orbs.mld`. Then specify this file in `dens_ana.in`:

::

    mo_file=orbs.mld

Firefly - TDDFT
~~~~~~~~~~~~~~~

Firefly has been succesfully interfaced with TheoDORE, see [EXAMPLES/SnH4-ecp.firefly](https://sourceforge.net/p/theodore-qc/code/ci/master/tree/EXAMPLES/SnH4-ecp.firefly/). Firefly output is parsed with the [cclib library](http://cclib.github.io/).

ADF - TDDFT
~~~~~~~~~~~
In the new ADF interface all information is read from the binary **TAPE21** file. Use the <code>rfile</code> option to point to this file.

To run the analysis, you need to activate the ADF scripts and license, e.g.

::

    . ~/adfrc.sh
    export SCMLICENSE=/usr/license/adf/licenses/license.txt
    export PYTHONPATH=$PYTHONPATH:/usr/license/adf/adf2016.101/scripting

The interface analyzes the eigenvectors of the reduced dimensional problem as printed out by ADF. Note, that these are only normalized in the case of the Tamm-Dancoff approximation.

It is not possible to visualize NTOs using the TheoDORE/ADF interface since Slater type orbitals, as employed by ADF, are not supported. It is, however, possible to compute NTOs within ADF itself.

The atom-numbering for <code>at_lists</code> pertains to the original ordering in the input file rather than the internal ordering used by ADF.

DFTB+ - TDDFTB
~~~~~~~~~~~~~~

An interface to DFTB+ was written by Ljiljana Stojanovic. This interface currently reads the following files:

*  EXC.DAT (main excited state information) - specified as 'rfile'
*  eigenvec.out (MO coefficients) - specified as 'mo_file'
*  XplusY.DAT (response vector)
*  SPX.DAT (ordering of response vector)
*  geom.xyz (geometry information)
*  detailed.out (orbital occupations and energies)
*  wfc.3ob-3-1.hsd (DFTB parameter file)

Other programs (cclib)
~~~~~~~~~~~~~~~~~~~~~~
In principle all third party programs, which are parsed by the [cclib library](http://cclib.github.io/) can be used. These are: ADF, Firefly, GAMESS, Gaussian, Jaguar, Molpro, ORCA. But not all of these have been tested by the developers and it may be necessary to set some additional program specific options. Please report, if you did so successfully.

To quickly check whether a logfile can be parsed by cclib, simply type:

::

    cc_check.py <logfile>
