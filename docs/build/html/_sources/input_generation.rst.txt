Input generation
----------------

The `theoinp` utility allows quick generation of the input information needed by TheoDORE. When running `theoinp`, the relevant options are suggested by default and only minimal input by the user is needed. It is however advised to take a look at the [program specific information](Program specific information) for the different interfaced quantum chemistry codes. After running `theoinp` a list of [keywords](Keywords) is written to the file `dens_ana.in`.



.. colt_commandline:: theodore  run
   :subparsers: options(theoinp)
   :header: False

   alias = theodore

   [arg_format]
   name = 25
   comment = 60

   [subparser_format]
   name = 25
   comment = 60


Molecular fragment definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A key input quantity for the charge transfer number analysis in TheoDORE is the molecular fragment definition. First, the user has to decide on how to separate the system under study into fragments. There is no unique way to do so, and it may be necessary to try out different fragmentation schemes to produce the most meaningful results.

There are five different ways of specifying the fragment information in `theoinp`.

::

    Mode for specifying molecular fragments (at_lists):
        [ 1] Manual input
        [ 2] Automatic generation by fragment (using python-openbabel)
        [ 3] Automatic generation for transition metal complexes (using python-openbabel)
        [ 4] Automatic generation by element (using python-openbabel)
        [ 5] Leave empty and fill out later

In mode 1 you will be asked to enter the indices of the atoms that belong to the different fragments successively.

Mode 2 is a special utility for automatic fragment definition. If the system of interest is composed of different molecules, these are detected automatically. Further customization can be performed by exporting the molecule in `.mol` format and changing the bond definitions by using for example [Avogadro](http://avogadro.cc/).

Mode 3 is a shortcut version of mode 2 that specifically works for transition metal complexes, see [*Coord. Chem. Rev.*, **361**, 74 (2018)](http://dx.doi.org/10.1016/j.ccr.2018.01.019). You simply have to add the atom of the transition metal and the system is automatically separated into the transition metal and the different ligands.

Mode 4 automatically separates the list by elements.

For mode 5 the `dens_ana.in` file has to be edited manually. For example, the input

::

    at_lists = [ [1,3,4], [2,5,6] ]

means that there are two molecular fragments. The first contains atoms 1, 3, and 4, the second 2, 5, and 6.
