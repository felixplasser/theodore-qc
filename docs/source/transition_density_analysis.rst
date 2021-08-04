Transition density matrix analysis
----------------------------------

.. colt_commandline:: theodore run
   :subparsers: options(analyze_tden)
   :header: False

    main_order = usage, pos_args, opt_args, subparser_args, space
    alias = theodore

    [arg_format]
    name = 20
    comment = 60



Analysis of the transition density matrix (1TDM) is performed with the script `analyze_tden.py`. The main types of analysis that can be performed are charge transfer numbers, natural transition orbitals, and an exciton size analysis.

**Note:** Traditionally, only spin-restricted computations are supported by TheoDORE and the treatment of unrestricted computations via the new script `analyze_tden_unr.py`  (see below) is not well-tested yet.

Charge transfer numbers
~~~~~~~~~~~~~~~~~~~~~~~

The charge transfer numbers are computed as partial summations over squared transition density matrix elements of molecular fragments. The key step for a charge transfer analysis is to [divde the molecule into meaningful fragments](Input).

The actual charge transfer number analysis is a generalized population analysis. Currently, three formulas are implemented: `Om_formula=0` yields a Mulliken style analysis

.. math::
   
    \begin{eqnarray}
    \Omega_{AB}=\sum_{\mu\in A}\sum_{\nu\in B}\left(\mathbf{D}^{0I}\mathbf{S}\right)_{\mu\nu}\left(\mathbf{S}\mathbf{D}^{0I}\right)_{\mu\nu}
    \end{eqnarray}

:math:`Om\_formula=1` yields a somewhat different Mulliken style analysis

.. math::

    \begin{eqnarray}
    \Omega_{AB}=0.5\sum_{\mu\in A}\sum_{\nu\in B}\left(\mathbf{D}^{0I}\mathbf{S}\right)_{\mu\nu}\left(\mathbf{S}\mathbf{D}^{0I}\right)_{\mu\nu}+D^{0I}_{\mu\nu}\left(\mathbf{S}\mathbf{D}^{0I}\mathbf{S}\right)_{\mu\nu}
    \end{eqnarray}

:math:`Om\_formula=2` yields a Lowdin style analysis

.. math::

    \begin{eqnarray}
    \Omega_{AB}=\sum_{\mu\in A}\sum_{\nu\in B}\left(\mathbf{S}^{1/2}\mathbf{D}^{0I}\mathbf{S}^{1/2}\right)_{\mu\nu}{}^{2}
    \end{eqnarray}

The last case (:math:`Om\_formula=2`) is probably the best option in terms of, both, computational effort and numerical stability.

Subsequently TheoDORE computes fragment based descriptors of the &Omega;-matrix. These are specified using the `prop_list` keyword, e.g.

::

    prop_list=['Om', 'POS', 'PR', 'CT', 'CTnt']

The complete list of available descriptors can be found [here](Transition density matrix analysis/attachment/Om_desc.pdf), see also the source file [Om_descriptors.py](https://sourceforge.net/p/theodore-qc/code/ci/master/tree/lib/Om_descriptors.py).

Natural transition orbitals
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The natural transition orbitals (NTOs) are constructed through a singular value decomposition of the transition density matrx **D**<sup>0I</sup>

.. math::

    \mathbf{D}^{0I}=\mathbf{U}diag(\lambda_1,\lambda_2,\ldots \mathbf{V}^T)


The are written to either a Molden format file or a compact script, which can be executed in Jmol. In the case of the Molden file the following convention is adopted: NTO singular values are written into the Occ field of the Molden file where negative values denote hole orbitals and positive values electron orbitals.

To count the number of NTO transitions involved, use the NTO participation ratio PR<sub>NTO</sub>

.. math::
    PR_{NTO}=\frac{\Omega}{\sum_i\lambda_i^2}

or the entanglement measures defined in [*J. Chem. Phys.*, **144**, 194107 (2016)](http://dx.doi.org/10.1063/1.4949535)

.. math::
    S_{H|E}=-\sum_i\lambda_i\log_2\lambda_i

and

.. math::
    \begin{eqnarray}
    Z_{HE}=2^{S_{H|E}},
    \end{eqnarray}

which are accessed by the following keywords:

::

    prop_list=['PRNTO', 'S_HE', 'Z_HE']

The *hole*/*particle* densities, which are analogous to the [attachment/detachment densities](State density matrix analysis), can be obtained as weighted sums over squared NTOs. This has to be [done externally](Orbitals and Densities), e.g. in Molden.

Conditional densities / domain NTOs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An alternative to the charge transfer numbers and NTOs has been introduced in [*ChemPhotoChem* (2019)](http://dx.doi.org/10.1002/cptc.201900014). To compute conditional densities or domain NTOs, specify in `theoinp`

::

    Perform analysis of domain NTOs and conditional densities? (comp_dntos):
    Choice (y/n): [n] y

Using [ORBKIT](https://orbkit.github.io/), you can directly generate cube files of the conditional densities, see also [Orbitals and Densities](Orbitals%20and%20Densities/). Specify whether you want to compute the conditional densities for fixed hole (1), fixed electron (2) or both (3).

::

    Compute conditional densities as cube files?
     0 - no, 1 - hole, 2 - electron, 3 - both (comp_dnto_dens):
    Choice: [0] 1

Plotting in VMD using [vmd_plots.py](Utility scripts/#vmd_plotspy).

Exciton size analysis
~~~~~~~~~~~~~~~~~~~~~

An approximate exciton size, [*PCCP*, **18**, 2548 (2016)](http://dx.doi.org/10.1039/c5cp07077e), (computed as the root-mean-square *electron-hole* separation, denoted `RMSeh`) is constructed as

.. math::
    d_{exc}=\sqrt{\sum_{MN}\Omega_{MN}d_{MN}^2/\Omega}

where M and N are two atom indices and d<sub>MN</sub> is the distance between them. The result is given in Angstrom.

::

    prop_list=['RMSeh']

Analysis of unrestricted computations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unrestricted computations are supported via the script `analyze_tden_unr.py`. This has only been tested with ORCA. For a well-tested support of unrestricted computations you have to resort to the implementations in Q-Chem and OpenMolcas.

The script `analyze_tden_unr.py` performs independent calculations for alpha and beta spin and writes the results to the subdirectories `ALPHA` and `BETA`. Natural transition orbitals can be written into these subdirectories as Molden files. Subsequently, the information is added up and collected in the main directory.

Analysis of spin-orbit coupled states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The analysis of spin-orbit coupled states, see [*Coord. Chem. Rev.*, **361**, 74 (2018)](http://dx.doi.org/10.1016/j.ccr.2018.01.019), is possible using the script `analyze_tden_soc.py`. Note, however, that this analysis is still in an experimental stage and is only possible for ADF.
