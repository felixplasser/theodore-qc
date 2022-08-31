Plotting
--------

Electron-hole correlation plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Electron-hole* correlation plots can be created with the ``theodore plot_omfrag`` utility.
You are asked about the output options directly by this script.
To view the output, simply open the file ``OmFrag.html`` in a web browser (`Example plots <http://theodore-qc.sourceforge.net/images/OmFrag.html>`_).

Generation of property graphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``plot_graph`` allows automatic generation of graphs of the different properties evaluated by TheoDORE.
For this purpose, it is assumed that several directories with analagous computations are present (e.g. from a potential curve).

Call ``theodore plot_graph`` and follow the instructions on the screen.

The graphs are collected in a file ``graphs.html`` (`Example graphs <http://theodore-qc.sourceforge.net/images/graphs.html>`_).

For parsing `Newton-X <http://www.newtonx.org/>`_ trajectories, a specialized script ``plot_graph_nx.py`` exists.

Absorption spectrum
~~~~~~~~~~~~~~~~~~~
Create a convoluted absorption spectrum from the energies and oscillator strengths parsed by ``analyze_tden``.

::

    theodore spectrum <tden_summ1> [<tden_summ2> ...]

It is also possible to compute the density of states (no weighting by oscillator strengths) and to add restrictions with respect to the states chosen (e.g. only states with CT > 0.5).

Fragment decomposition (e/h populations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An interactive script ``plot_frag_decomp`` can be used for plotting the fragment decomposition of the electron and hole populations of the excited states.
To be called after ``analyze_tden``

.. figure:: figures/frag_decomp.png

Fragment decomposition (Omega matrices)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An interactive script ``plot_om_bars`` can be used for plotting the fragment decomposition of the electron and hole populations of the excited states.
To be called after ``analyze_tden``.
The script proceeds by asking for different components, which are each in turn composed of several donor/acceptor pairs.
To define, for example, the MLCT states in a complex with three ligands, proceed like this:

.. code-block:: text

    Name of component 1
    Choice: MLCT

    Color for plotting
    Choice: blue

    *** Fragment pairs belonging to MLCT ***
     Enter two indices between 1 and 4, separated by spaces
     Leave empty to finish

    Hole/electron indices for pair 1
    Choice: 1 2

    Hole/electron indices for pair 2
    Choice: 1 3

    Hole/electron indices for pair 3
    Choice: 1 4

    Hole/electron indices for pair 4
    Choice: 

    ... switching to next component.

.. figure:: figures/Om_bars.png
