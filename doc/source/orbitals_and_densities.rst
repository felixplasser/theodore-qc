Orbitals and Densities
----------------------


Orbitals are exported in three forms:

* [Molden format](http://www.cmbi.ru.nl/molden/molden_format.html), which can be directly read by a number of visualization programs

* As a compact script for the [Jmol](http://jmol.sourceforge.net/) program

* Cube files of densities created by [ORBKIT](http://orbkit.github.io/)

Using the Jmol script
~~~~~~~~~~~~~~~~~~~~~

After running <code>analyze_tden.py</code> (<code>analyze_sden.py</code>) the file <code>nto_jmol.spt</code> (<code>ndo_jmol.spt</code>) is created. This is an input file for jmol, which can be directly executed as

::

    jmol -n nto_jmol.spt

For more flexibility, the following semi-interactive procedure can be applied:

1. Open Jmol

+ Run the first few lines of the script to preview the settings (simply copy them into the Jmol console)

+ Adjust the perspective

+ Run the remaining lines of the script (by copying into the Jmol console)

+ Open the file <code>nto.html</code> (<code>ndo.html</code>) to view the result ([Example](http://theodore-qc.sourceforge.net/images/nto.html)).

Density plotting (ORBKIT)
~~~~~~~~~~~~~~~~~~~~~~~~~

Cube files of densities can be directly created with [ORBKIT](http://orbkit.github.io/). The interface is controlled via <code>theoinp</code>.

The cube files can be loaded into VMD and visualized using VMD network files provided by ORBKIT or by using the [vmd_plots.py](Utility%20scripts/#vmd_plotspy) facility of TheoDORE. For the latter, just run:

::

    vmd_plots.py *.cb

ORBKIT is somewhat sensitive in terms of the molden file with orbital information. To achieve the best result, it is advisable to create this file within Molden rather than through TheoDORE.

Density plotting (Molden)
~~~~~~~~~~~~~~~~~~~~~~~~~

The weighted sum over the NTOs gives the *hole* / *particle* densities, while the NDOs yield the *attachment* / *detachment* densities.

The following procedure is suggested:

* Create Molden format files of the NTOs / NDOs using analyze_tden.py / analyze_sden.py

+ To extract the separate *hole* / *particle* contributions use the [extract_molden.py script](Utility scripts)

Run, e.g.:

::

    extract_molden.py nto_1-1-a.mld

This will create a new directory <code>nto_1-1-a.mld.dir</code> containing the files <code>nto_1-1-a.mld_elec.mld</code> and
<code>nto_1-1-a.mld_hole.mld</code> with the separate contributions.

* Open <code>nto_1-1-a.mld_elec.mld</code> in Molden (or preferably gmolden)
    * Go to "density mode"
    * Select "Plot Function" - "Density"
    * Use "Plot Mode" - "Space"
    * The result is the *particle* density as defined in [this Ref](http://dx.doi.org/10.1063/1.4885819).
