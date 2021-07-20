Installation
------------

The installation on a Linux system is split into four simple parts:

Download and extract the source file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the newest release TheoDORE_v.s.tgz from the [download page](/p/theodore-qc/files/).

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

Starting with TheoDORE 2.0, the suggested python version is python 3.5 but TheoDORE is also compatible with python 2.7.14 .

*Note*: Due to the switch to python3 it is necessary that all external packages are available under python3.

External packages
~~~~~~~~~~~~~~~~~

Four external packages are used by TheoDORE:


    - [python3-numpy](http://numpy.scipy.org/) - for basic numerical manipulations
    - [python3-matplotlib](http://matplotlib.sourceforge.net/) (optional) - for plotting of graphs
    - [python-openbabel](http://openbabel.org/wiki/Python) (optional) - for extended file-parsing capabilities of molecular structure files
    - [ORBKIT](http://orbkit.github.io/) (optional) - for creating cube files of densities

The `numpy` and `matplotlib` packages are usually readily available with the standard installation tools, e.g. `apt-get`, `yum` etc. Alternatively, they may be downloaded from the URLs specified. If no integrated installation is performed, then it is necessary to add these libraries to the `PYTHONPATH` (see above).

For the installation of `orbkit` use the dedicated [TheoDORE fork on github](https://github.com/felixplasser/orbkit) and follow the installation instructions given there.

*Note*: the [cclib](http://cclib.github.io/) package is distributed as part of TheoDORE (starting in Version 1.3) and no separate installation step is necessary.

Using anaconda
~~~~~~~~~~~~~~

A straightforward and universal way of installing most of the required packages is through the use of Anaconda.

First download the [anaconda distribution](https://www.anaconda.com/distribution/) and do the installation. Then run the commands:

::

    conda install numpy matplotlib
    conda install -c openbabel openbabel

Testing
~~~~~~~

Bash script
___________

Use the `theo_test.bash` utility to test you installation:

::

    theo_test.bash [<module>]
      available modules: standard, all, openbabel, cclib, adf, noadf

Ideally the test routines should finish with the following lines:

::

 *** All tests finished (number of errors: 0)


pytest
______

The tests can also be invoked with pytest. In the `EXAMPLES` directory run

::

    pytest -v && cat pytest.out

If a test fails, run `pytest -s` to get more detailed information.

Installation on Windows 7
~~~~~~~~~~~~~~~~~~~~~~~~~
(by Dr. Siddheshwar Chopra)

1. First ensure that you have Python 2.7.9 or above installed.

2. Download the required packages (ending with .whl) for your OS (32 or 64 bit) and Python version (e.g. cp27 for Python 2.7) from the website of Christoph Gohlke [http://www.lfd.uci.edu/~gohlke/pythonlibs/](http://www.lfd.uci.edu/~gohlke/pythonlibs/). Download at least NumPy, SciPy, Matplotlib.
3. Now install the packages (\*.whl) by typing "pip install packagename" without quotes.
4. "cclib" also could be installed. You can get cclib-1.3 for the python version. And follow the installation steps.
5. You might also need to install python-openbabel. (openbabel-python-1.6.py27.exe)
6. Set the environment variables:

::

    THEODIR= C:\TheoDORE_1.1.1 (TheoDORE_1.1.1 installation folder path)
    PATH=C:\TheoDORE_1.1.1\bin      (#binary files path)
    PYTHONPATH=C:\TheoDORE_1.1.1\lib  (#library files path)  

Note: In Windows 7, you can add/edit the environment variables as follows:

a) Right click on "My Computer". Select "Properties".
b) Select "Advanced System Settings" from the left column.
c) Select "Environment Variables" from "Advanced" Tab.
d) Here you can click "Add" or "Edit" to add these variables. You might need to restart.

Ensure to make the installation folder of TheoDORE, without any spaces in its name. Example "C:/TheoDORE_1.1.1" is appropriate.
