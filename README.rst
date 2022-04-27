TheoDORE
--------

The **TheoDORE** (Theoretical Density, Orbital Relaxation and Exciton analysis) package is a program suite for the analysis of excited states obtained from quantum chemical excited state calculations.

*Author*: Felix Plasser

*Contributors*: Ljiljana Stojanovic, Gunter Hermann, Sebastian Mai, Maximilian F.S.J. Menger, Patrick Kimber

TheoDORE is distributed under the GNU General Public License 3.0 (see `LICENSE.txt <https://github.com/felixplasser/theodore-qc/blob/master/LICENSE.txt>`_).

.. image:: ./doc/source/_static/theodore.png

Documentation
~~~~~~~~~~~~~
* `Documentation of TheoDORE 2 <https://sourceforge.net/p/theodore-qc/wiki/Home/>`_
* `Documentation of TheoDORE 3 <https://theodore-qc.sourceforge.io/doc_theo3-alpha/index.html>`_

For bugs / feature requests, please use the `issues page <https://github.com/felixplasser/theodore-qc/issues>`_.

Installation
~~~~~~~~~~~~
* You can obtain the newest TheoDORE release from `github releases <https://github.com/felixplasser/theodore-qc/releases>`_.
* Alternatively, to get the current development version of the code and test suite, run

::

    git clone --recursive https://github.com/felixplasser/theodore-qc.git
    git clone https://github.com/felixplasser/theo-test.git

* To run TheoDORE, setup ``PATH`` and ``PYTHONPATH`` as explained in the `installation instructions <https://theodore-qc.sourceforge.io/doc_theo3-alpha/installation.html>`_.

Usage
~~~~~
A detailed description of the usage is given `here <https://theodore-qc.sourceforge.io/doc_theo3-alpha/usage.html>`_.

* As opposed to earlier versions, TheoDORE 3 is activated with a central driver script ``theodore``.
* To get a list of all implemented options simply type

::

    theodore -h

* For input generation

::

    theodore theoinp

* Analysis of transition density matrices

::

    theodore analyze_tden

External libraries
~~~~~~~~~~~~~~~~~~

TheoDORE uses the cclib library (http://cclib.github.io) for some of its file parsing work.
cclib is distributed under a BSD 3-Clause License, see cclib/LICENSE .
Copyright (c) 2017, the cclib development team.
cclib is contained as part of this distribution.

TheoDORE uses colt (https://github.com/mfsjmenger/colt) for its commandline interface.
colt is distributed under the Apache License 2.0.
colt is contained as part of this distribution.

Citation for cclib:

``N. M. O'Boyle, A. L. Tenderholt, K. M. Langner, J. Comput. Chem. (2008), 29, 839.``

TheoDORE uses Open Babel (http://openbabel.org/) for reading chemical structure files.
Open Babel is distributed under a GPL license. Install Open Babel if you want to use this functionality.
Citations for Open Babel:

``N. M. O'Boyle, M. Banck, C. A. James, C. Morley, T. Vandermeersch, and G. R. Hutchison, J. Cheminf. (2011), 3, 33.``

``The Open Babel Package, version 2.3.1 http://openbabel.org``

TheoDORE uses ORBKIT (http://orbkit.github.io/) for visualization of orbitals and densities.
ORBKIT is distributed under an LGPL license. Install orbkit if you want to use this functionality.
Citation for ORBKIT:

``G. Hermann, V. Pohl, J. C. Tremblay, B. Paulus, H.-C. Hege, A. Schild, J. Comput. Chem. (2016), 37, 1511.``

Disclaimer
~~~~~~~~~~

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
