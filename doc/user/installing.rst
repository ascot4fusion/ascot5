.. _Installing:

============
Installation
============

.. tab-set::

   .. tab-item:: Full installation

      .. card::

         The most convenient way to install ASCOT5 is to use `Conda <https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html>`_.

         **Without MPI:**

         .. code-block:: bash

            git clone https://github.com/ascot4fusion/ascot5.git
            cd ascot5
            conda env create -f environment.yaml
            conda activate ascot-env
            make ascot5
            pip install -e .

         The executables are located in the ``build`` directory.

         **With MPI:**

         .. code-block:: bash

            git clone https://github.com/ascot4fusion/ascot5.git
            cd ascot5
            conda env create -f environment-mpi.yaml
            conda activate ascot-mpi
            make libascot -j MPI=1
            make ascot5_main -j MPI=1 # Also bbnbi5 if needed
            pip install -e .

         Do note that this method uses the MPI packaged with Conda.
         This might not be preferred on some HPC clusters, so it is usually best to use the native MPI library instead.
         `To switch to native MPI <https://conda-forge.org/docs/user/tipsandtricks/#using-external-message-passing-interface-mpi-libraries>`_, do the following:

         .. code-block:: bash

            module load <local-mpi-library>
            mpirun -V # Prints the version of the MPI library to be used in the next line
            conda install "<local-mpi-library>=x.y.z=external_*"
            export HDF5_CC=$(which mpicc)
            export HDF5_CLINKER=$(which mpicc)

            # If you have installed mpi4py already, it needs to be reinstalled
            # with the new MPI library.
            python -m pip cache remove mpi4py
            python -m pip install mpi4py

         Conda does not have all versions of MPI libraries available, so one might have to use asterix in the minor version number, for example:

         .. code-block:: bash

            conda install "openmpi=4.*=external_*"

         .. note::
            Add the following lines in your batch script to use the conda environment when submitting jobs via SLURM:

            .. code-block:: bash

               eval "$(conda shell.bash hook)"
               conda activate ascot-env

            If Conda is not available on your cluster, it can be easily `installed <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_ (doesn't require sudo) with:

            .. code-block:: bash

               curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
               bash Miniforge3-$(uname)-$(uname -m).sh


   .. tab-item:: Minimal installation (for HPC usage)

      .. card::

         Building the code from the source (without Conda) provides a light installation with the capability to run simulations using native libraries and perform limited pre- and post-processing.
         A typical way to use ASCOT5 is to run simulations on a computing cluster and carry out pre- and postprocessing on a home cluster or a workstation.
         For optimal performance, use this method on the HPC cluster for running simulations, and then the full installation where you process the data.

         .. rubric:: Requirements

         1. Install the requirements or use the module system:

            - C compiler (Intel)
            - HDF5
            - OpenMP
            - MPI
            - Python3.10

         2. Download the source and set up the virtual environment:

            .. code-block:: bash

               git clone https://github.com/ascot4fusion/ascot5.git
               python -m venv ascot-env
               source activate ascot-env/bin/activate

         3. Install ``a5py`` and compile the executables which will be located at ``build/``:

            .. code-block:: bash

               cd ascot5
               pip install -e .
               make ascot5_main -j MPI=1
               make libascot -j MPI=1

         See :ref:`here<Compiling>` for tips on how to compile the code on different platforms.

         *Optional* Add the following lines to your `.bashrc` to automatically activate ASCOT5 environment each time you login:

         .. code-block::

            <module loads and exports here>
            source activate /path/to/ascot5env

         *GPU* To compile the code for the NVIDIA GPU nodes, you'll need ``nvc`` compiler:

         .. code-block:: bash

            make ascot5_main -j GPU=1 ACC=1 CC=nvc

         Currently AMD GPUs are not supported.

   .. tab-item:: Developers

      .. card::

         This is a full installation from the source (using Conda) and with the optional packages present that are required to build ``ascot2py.py`` and the documentation.
         Note that the first step requires you to `add SSH keys on GitHub <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_ whenever on a new machine.

         .. code-block:: bash

            git clone git@github.com:ascot4fusion/ascot5.git
            cd ascot5
            conda env create -f environment-dev.yaml
            conda activate ascot-dev
            make libascot -j
            make ascot5_main -j
            pip install -e .

The simulation options are edited with the local text editor, which usually happens to be ``vim``.
Consider adding the following line to your ``.bashrc`` (or ``.bash_profile`` if working locally):

.. code-block:: bash

   export EDITOR=/usr/bin/emacs

Whenever there is a new release, you can update the code as:

.. code-block:: bash

   git pull
   make clean
   make ascot5_main -j (MPI=1)
   make libascot -j (MPI=1)
   pip install -e .

Always use the ``main`` branch when running simulations unless you specifically need something from a feature branch.
Version numbers are specified via tags.
To switch to a different version number:

.. code-block:: bash

   git checkout <version-number> # e.g. 5.5.3
   make clean
   make ascot5_main -j (MPI=1)
   make libascot -j (MPI=1)
   pip install -e .

Test that ASCOT5 was properly installed by running the :ref:`introduction<Tutorial>`.

.. admonition:: Troubleshooting

   **Conda complains that** ``freeqdsk`` **was not found or could not be installed.**

   - Remove ``freeqdsk`` from ``environment.yaml`` and install it with ``pip`` once the environment is activated.
     The cause of this error is unknown.

   **Make stops since** ``libhdf5*.a`` **could not be located.**

   - The compiler should be using shared libraries, not static.
     Add the following flag for the compiler: ``make ascot5_main FLAGS="-shlib"``

.. _Compiling:

Compiling on different platforms
================================

See host-specifics

.. _Compilerflags:

Settings when compiling
=======================

Some of the ASCOT5 options require recompiling the code.
Parameters that can be given arguments for ``make`` are (the default values are shown)

.. code-block:: bash

   make -j ascot5_main NSIMD=16 CC=h5cc

.. list-table::
   :widths: 10 50

   * - NSIMD
     - Number of particles simulated in parallel in each SIMD vector.
       These are processed simultaneously by each thread and the optimal number depends on the hardware.
       If unsure, keep the default value.
   * - CC
     - C compiler to use.

Compiler flags can be provided with ``FLAGS`` (and linker flags with ``LFLAGS``) parameter, e.g.

.. code-block:: bash

   make -j ascot5_main FLAGS="-qno-offload"

Some parameters relevant for ASCOT5 are (these are compiler dependent):

.. list-table::
   :widths: 10 50

   * - ``-qno-openmp-offload`` or ``-foffload=disable``
     - Disables offload.
       Recommended when not using Xeon Phi.
   * - ``-diag-disable 3180``
     - Disables Intel compiler warnings about unrecognized pragmas when the offloading is disabled.
   * - ``-xcommon-avx512``, ``-xcore-avx512``, ``-xmic-avx512``
     - Compile the code for Skylake or KNL processors, optimize for Skylake, optimize for KNL.
   * - ``-vecabi=cmdtarget``
     - Enables vector instructions for NSIMD > 2.
   * - ``-ipo``
     - "Interprocedural Optimization" which might increase the performance somewhat.
   * - ``-qopt-report=5`` and ``-qopt-report-phase=vec``
     - Generate vectorization reports in \*optrpt files.
       Only useful for developers.

Additional compile-time parameters can be found in ``ascot5.h``, but there is rarely a need to change those.


For Developers
==============

There is a separate Conda environment for setting up virtual environment.


To build user documentation

or developer documentation

or both

and these are located in


Tests are run with

Tests are implemented using pytest so you can also run individual tests with


Typehints are checked with mypy

or

to check individual files.

Linting is done with pylint

or

for individual files.
