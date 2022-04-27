.. _orb-dens:
Orbitals and Densities
----------------------

Orbitals are exported in three forms:

* `Molden format <http://www.cmbi.ru.nl/molden/molden_format.html>`_, which can be directly read by a number of visualization programs

* As a compact script for the `Jmol <http://jmol.sourceforge.net/>`_ program

* Cube files of densities created by the `ORBKIT <http://orbkit.github.io/>`_

Using the Jmol script
~~~~~~~~~~~~~~~~~~~~~

After running ``analyze_tden`` (``analyze_sden``) the file ``nto_jmol.spt`` (``ndo_jmol.spt``) is created.
This is an input file for jmol, which can be directly executed as

::

    jmol -n nto_jmol.spt

For more flexibility, the following semi-interactive procedure can be applied:

+ Open Jmol

+ Run the first few lines of the script to preview the settings (simply copy them into the Jmol console)

+ Adjust the perspective

+ Run the remaining lines of the script (by copying into the Jmol console)

+ Open the file ``nto.html`` (``ndo.html``) to view the result (`Example <http://theodore-qc.sourceforge.net/images/nto.html>`_).

Density plotting (ORBKIT)
~~~~~~~~~~~~~~~~~~~~~~~~~

Cube files of densities can be directly created with `ORBKIT <http://orbkit.github.io/>`_. The interface is controlled via ``theoinp``.

The cube files can be loaded into VMD and visualized using VMD network files provided by ORBKIT or by using the ``vmd_plots`` facility of TheoDORE. For the latter, just run:

::

    theodore vmd_plots *.cb

ORBKIT is somewhat sensitive in terms of the molden file with orbital information. To achieve the best result, it is advisable to create this file within Molden rather than through TheoDORE.

Density plotting (Jmol / Molden)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The weighted sum over the NTOs gives the *hole* / *particle* densities, while the NDOs yield the *attachment* / *detachment* densities.
First, create Molden format files of the NTOs / NDOs using analyze_tden.py / analyze_sden.py

Jmol
^^^^

To use Jmol use the ``jmol_mos`` tool. Run, e.g.

::

    theodore jmol_mos nto*.mld

Activate occupancy weighted densities using the following set of options:

::

    Specification of the orbital indices to be plotted (spec):
      [ 1]       sten - Start and end indices
      [ 2]   frontier - Number of frontier orbitals
      [ 3]        occ - Occupation threshold
    Choice: 3

    ...

    Preprocess and merge the Molden files (preprocess):
    Choice (y/n): [n] n

    Compute occupancy-weighted densities (attachment/dentachment, hole/electron, etc)? (do_dens):
    Choice (y/n): [n] y

This will plot densities corresponding to all orbitals that fit the occupation threshold.
Separate densities (``dens_plus.png``, ``dens_minus.png``) are computed for orbitals with positive and negative eigenvalues.

Molden
^^^^^^

For Molden the following procedure is suggested.

+ To extract the separate *hole* / *particle* contributions use the ``extract_molden`` functitonality

Run, e.g.:

::

    theodore extract_molden nto_1-1-a.mld

This will create a new directory ``nto_1-1-a.mld.dir`` containing the files ``nto_1-1-a.mld_elec.mld`` and
``nto_1-1-a.mld_hole.mld`` with the separate contributions.

* Open ``nto_1-1-a.mld_elec.mld`` in Molden (or preferably gmolden)
* Go to "density mode"
* Select "Plot Function" - "Density"
* Use "Plot Mode" - "Space"
* The result is the *particle* density as defined in `JCP, 141, 024106 (2014) <http://dx.doi.org/10.1063/1.4885819>`_.
