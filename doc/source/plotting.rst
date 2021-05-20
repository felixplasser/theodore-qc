Plotting
--------

Electron-hole correlation plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Electron-hole* correlation plots can be created with the <code>plot_OmFrag.py</code> utility. You are asked about the output options directly by this script. To view the output, simply open the file <code>OmFrag.html</code> in a web browser ([Example](http://theodore-qc.sourceforge.net/images/OmFrag.html)).

Generation of property graphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The <code>plot_graph.py</code> allows automatic generation of graphs of the different properties evaluated by TheoDORE. For this purpose, it is assumed that several directories with analagous computations are present (e.g. from a potential curve).

Call <code>plot_graph.py</code> and follow the instructions on the screen.

The graphs are collected in a file <code>graphs.html</code> ([Example](http://theodore-qc.sourceforge.net/images/graphs.html)).

For parsing [Newton-X](http://www.newtonx.org/) trajectories, a specialized script <code>plot_graph_nx.py</code> exists.

Absorption spectrum
~~~~~~~~~~~~~~~~~~~
Create a convoluted absorption spectrum from the energies and oscillator strengths parsed by <code>analyze_tden.py</code>.

::

    spectrum.py <tden_summ1> [<tden_summ2> ...]

It is also possible to compute the density of states (no weighting by oscillator strengths) and to add restrictions with respect to the states chosen (e.g. only states with CT > 0.5).

Fragment decomposition (e/h populations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An interactive script <code>plot_frag_decomp.py</code> can be used for plotting the fragment decomposition of the electron and hole populations of the excited states. To be called after <code>analyze_tden.py</code>

![frag_decomp.png](https://sourceforge.net/p/theodore-qc/wiki/Plotting/attachment/frag_decomp.png)

Fragment decomposition (Omega matrices)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An interactive script <code>plot_Om_bars.py</code> can be used for plotting the fragment decomposition of the electron and hole populations of the excited states. To be called after <code>analyze_tden.py</code>. The script proceeds by asking for different components, which are each in turn composed of several donor/acceptor pairs. To define, for example, the MLCT states in a complex with three ligands, proceed like this:

::

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

![Om_bars.png](https://sourceforge.net/p/theodore-qc/wiki/Plotting/attachment/Om_bars.png)