State density matrix analysis
-----------------------------

Analysis of state and difference density matrices is performed with the script <code>analyze_sden.py</code>. In the case of an analysis of natural orbital (NO) files, also the <code>analyze_NOs.py</code> may be applied.

Population analysis
~~~~~~~~~~~~~~~~~~~

TheoDORE features Mulliken population analysis capabilities, which can be applied to normal densities, densities of effectively unpaired electrons or attachment/detachment densities (see below).

Unpaired electrons
~~~~~~~~~~~~~~~~~~

Two measures for a number of effectively unpaired electrons can be computed:

.. math::
    n_u=\sum_i\min(n_i,2-n_i)

and

.. math::
    n_{u,nl}=\sum_i n_i^2(2-n_i^2)

where *n<sub>i</sub>* marks the occupation number of natural orbital *i*.

Attachment/detachment analysis and natural difference orbitals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The difference density between an excited and a reference state can be diagonalized to provide the natural difference orbitals (NDOs). Summing over all NDOs of positive (negative) sign, the attachment (detachment) densities can be constructed. This has to be performed by an [external program](Orbitals and Densities).

Bond order and valence analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bond order and valences are computed according to the formulas given by I. Mayer, see [Int. J. Quant. Chem. **1986** XXIX, 477](http://dx.doi.org/10.1002/qua.560290320) for more information.