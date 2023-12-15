ASCOT5
======

`Documentation <https://ascot4fusion.github.io/ascot5/>`_

High-performance orbit-following code for fusion plasma physics and engineering.

New (and old) users are welcome to join our weekly meetings or Slack channel to discuss their research and plans with regards to ASCOT5.

This repository is maintained by ASCOT team in Aalto University and VTT Technical Research Centre of Finland.

Installation
============

Clone the repository:

.. code-block:: bash

   git clone git@github.com:ascot4fusion/ascot5.git

.. rubric:: Requirements

- C compiler
- HDF5
- OpenMP
- Python >= 3.10 (pre- and postprocessing)
- MPI (optional)
- VTK (optional, for 3D wall visualization)

Minimal Installation
********************

For running ASCOT5 on this platform and performing pre- and post-processing on another platform:

.. code-block:: bash

   cd ascot5
   make ascot5_main

The binary is located at ``build/ascot5_main``.

Full Installation
*****************

For full installation both ascot5_main and libascot.so are needed:

.. code-block:: bash

   cd ascot5
   make ascot5_main
   make libascot

Python Library (a5py)
*********************

Useful even for minimal installation; it provides command-line tools for updating simulation options.

Create a virtual environment (optional but recommended), activate it, and install ``a5py``:

.. code-block:: bash

   cd ..
   virtualenv -p python3 --system-site-packages ascotenv
   source ascotenv/bin/activate
   pip install -e ascot5/

How to Contribute
=================

.. admonition:: As an User:

   - Verify your results and report violations in issues.
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
