ASCOT5
======

`Documentation <https://ascot4fusion.github.io/ascot5/>`_

High-performance orbit-following code for fusion plasma physics and engineering.

New (and old) users are welcome to join our weekly meetings or Slack channel to discuss their research and plans with regards to ASCOT5.

This repository is maintained by ASCOT team in Aalto University and VTT Technical Research Centre of Finland.

Installation
============

Most convenient way to install ASCOT5 is to use the Conda environment that comes with the source code.
Compile the main program and the library, and install the associated Python package with ``pip``.
Note that this installation is suitable mostly for pre- and postprocessing since it does not use MPI.
See the `documentation <https://ascot4fusion.github.io/ascot5/installing.html>`_ for more detailed instructions.

.. code-block:: bash

   git clone https://github.com/ascot4fusion/ascot5.git
   cd ascot5
   conda env create -f environment.yaml
   conda activate ascot-env
   make libascot -j
   make ascot5_main -j
   pip install -e .

How to Contribute
=================

.. admonition:: As an User:

   - Verify your results and report problems or physics violations in issues.
   - Add `compiling instructions <https://ascot4fusion.github.io/ascot5/installing.html#compiling-on-different-platforms>`_ for popular platforms and consider updating `the table with example simulation times <https://ascot4fusion.github.io/ascot5/simulations.html#examples>`_ for novel entries.
   - When `benchmarking <https://ascot4fusion.github.io/ascot5/testing.html#benchmarks>`_ ASCOT5 against other codes or `validating <https://ascot4fusion.github.io/ascot5/testing.html#validation>`_ it against experiments, please contact the maintainers to archive the simulation for use as a regression test.

.. admonition:: As a Developer:

   - Don't let the code daunt you!
     We're here to assist with any feature contributions, whether it's a small post-processing tool, a new plot, or an enhancement of an existing feature.
   - Start by creating an issue, then (fork and) make a branch ``feature/<issuenumber>-issue`` from develop.
   - When ready to merge, create a pull request running automated tests on your branch.
   - Upon test completion and acceptance, your feature will be merged into develop for inclusion in the next release.

Licence
=======

The ASCOT5 and associated programs are distributed under the terms of the GNU Lesser General Public License (LGPL).
Please see the files COPYING and COPYING.LESSER for more information.

This has been done after the code was released to the original authors by the Dean of School of Science of Aalto University and discussion between the key contributors, including Jari Varje, Konsta Särkimäki, Antti Snicker and Simppa Äkäslompolo.
