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
- Python3 (pre- and postprocessing)
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

GIT Usage
=========

When working on issues:

- Create a branch named ``feature/<issuenumber>-branchname`` for feature requests
- Create a branch named ``bugfix/<issuenumber>-branchname`` for bug fixes.
