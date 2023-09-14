ASCOT5
======

High-performance orbit-following code for fusion plasma physics and engineering.

If you are a new user, there is no need to ask for a permission to use the code.
However, we, the developers, would welcome it if you'd like to introduce yourself!
You are welcome to present the ASCOT5 work you have done or plan to do in one of our weekly meetings.

This repository is maintained by ASCOT team in Aalto University and VTT Technical Research Centre of Finland.

Installation
============

Minimal
*******

If you plan to run ASCOT5 on this platform and perform pre- and postprocessing on another platform, you can then make the minimum installation here.

.. code-block:: bash

   git clone git@github.com:ascot4fusion/ascot5.git
   cd ascot5
   make ascot5_main

The binary is located at ``build/ascot5_main``.

Python library a5py
*******************

The Python library is useful even if you do only the minimal installation as it provides command line tools to e.g. update simulation options.

It is recommended to set a virtual environment when working with ASCOT5.
You can use either ``virtualenv`` (3rd party tool)

.. code-block:: bash

   virtualenv -p python3 --system-site-packages ascotenv
   source ascotenv/bin/activate

or ``venv`` (packs with Python 3)

.. code-block:: bash

   python3 -m venv --system-site-packages ascotenv
   source ascotenv/bin/activate

or Conda (preferred)

.. code-block:: bash

   conda create --file ascot5/requirements.txt --name ascotenv
   conda activate ascotenv

Whichever was chosen, the python package can be installed with

.. code-block:: bash

   pip install -e ascot5

Here option ``-e`` makes any changes in the repository to be automatically reflected in your installation.

Full
****

In addition to above, the full installation requires library ``libascot.so`` since many of the features in ``a5py`` use the ASCOT5 C-routines directly.
The full installation enables simulations via Python and interpolating data exactly as it is interpolated in run time.

In ascot5 folder:

.. code-block:: bash

   make libascot


Documentation (for developers)
******************************

.. code-block:: bash

   cd ascot5
   make doc

The main page of the documentation is located at ``build/doc/index.html``.

Requirements
============

**C compiler**, **HDF5**, **OpenMP**, **Python3**, **MPI** (optional), and **VTK** (optional).
