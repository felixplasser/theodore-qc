Installation
------------

The installation on a Linux system is split into four simple parts:

Download and extract the source file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the newest release TheoDORE_v.s.tgz from the `github releases page <https://github.com/felixplasser/theodore-qc/releases>`_.

Extract the source file

::

    tar -xf TheoDORE_v.s.tgz

Setup the path specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To use TheoDORE the `PATH` and `PYTHONPATH` variables have to be adjusted. For this purpose you can use the provided bash script `setpaths.bash`

::

    #!/bin/bash
    export THEODIR=/yourpath/TheoDORE/TheoDORE_v.s
    export PATH=$THEODIR/bin:$PATH
    export PYTHONPATH=$THEODIR:$PYTHONPATH


**Note:** `PYTHONPATH` no longer points to the `lib` directory starting in TheoDORE_2.0.

Copy the above lines into your .bashrc file or run:

::

    source setpaths.bash

Alternatively, a csh script setpaths.csh is provided.

Python3
~~~~~~~

TheoDORE 3 is compatible with python3 (ideally >3.7).
No compatibility to python2 is maintained.

The older release, TheoDORE 2.4, is still compatible with python v2.7.14.

External packages
~~~~~~~~~~~~~~~~~

The following external packages are used by TheoDORE and require a separate installation:

    - `python3-numpy <http://numpy.scipy.org/>`_ - for basic numerical manipulations
    - `python3-matplotlib <http://matplotlib.sourceforge.net/>`_ *(optional)* - for plotting of graphs
    - `python3-openbabel <http://openbabel.org/wiki/Python>`_ *(optional)* - for extended file-parsing capabilities of molecular structure files
    - `ORBKIT <https://github.com/felixplasser/orbkit>`_ *(optional)* - For creating cube files of densities

    The first three are usually readily available with the standard installation tools, e.g. ``apt-get``, ``yum`` etc.
    Alternatively, they may be downloaded from the URLs specified.
    If no integrated installation is performed, then it is necessary to add these libraries to the `PYTHONPATH` (see above).
    For ORBKIT, please follow the `installation instructions <http://orbkit.github.io/install.html>`_.

The following external packages are provided along with the TheoDORE distribution and do not require installation.

    - `cclib <http://cclib.github.io/>`_ - for file parsing work
    - `periodictable <https://github.com/pkienzle/periodictable>`_ - information about periodic table (needed by cclib)
    - `colt <https://github.com/mfsjmenger/colt>`_ - user interface

Using anaconda
~~~~~~~~~~~~~~

A straightforward and universal way of installing most of the required packages is through the use of Anaconda.

First download the `anaconda distribution <https://www.anaconda.com/distribution/>`_ and do the installation. Then run the commands:

::

    conda install numpy matplotlib
    conda install -c openbabel openbabel

Testing
~~~~~~~

The tests are invoked with ``pytest-3``. In the ``EXAMPLES`` directory run (depending on the packages you want to test)

::

    pytest-3 UTILS STANDARD CCLIB
    pytest-3 EXTRA

The first three tests should always work.
For all tests in ``EXTRA`` to pass openbabel, ADF and ORBKIT need to be installed.

The tests are distributed within the standard TheoDORE distribution.
They are maintained via the `theodore-test <https://github.com/felixplasser/theodore-test>`_ repository.
