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
            make libascot -j
            make ascot5_main -j # Also bbnbi5 if needed
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

Here are instructions on how to compile ASCOT5 in some of the platforms where the code has been used.
If you are using ASCOT5 in a cluster that is not listed here, feel free to amend the list.

.. tab-set::

   .. tab-item:: CSC.fi puhti

      .. card::

         .. code-block:: bash

            module load StdEnv intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4 hdf5/1.10.4-mpi python-data

            make -j ascot5_main MPI=1

         Alternatively:

         .. code-block:: bash

            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"

   .. tab-item:: CSD3

      .. code-block:: bash

         module load hdf5/mpi/intel/2019.3/1.10.5
         make ascot5_main MPI=1 CC=mpicc FLAGS="-no-multibyte-chars -qno-openmp-offload -diag-disable 3180 -xmic-avx512 -vecabi=cmdtarget" LFLAGS="-lhdf5_hl -lhdf5" -j
         make libascot MPI=1 CC=mpicc FLAGS="-no-multibyte-chars -qno-openmp-offload -diag-disable 3180" LFLAGS="-lhdf5_hl -lhdf5" -j

   .. tab-item:: Freia (UKAEA)

      .. card::

         .. code-block:: bash

            make ascot5_main MPI=0 CC=gcc

         For libascot, one user needed to revert to python/3.5.1 and command

         .. code-block:: bash

            make libascot MPI=0

   .. tab-item:: ITER sdcc

      .. card::

	 For parallel version with MPI:

         .. code-block:: bash

            module load zlib/1.3.1-GCCcore-13.3.0 Szip/2.1.1-GCCcore-13.2.0  HDF5/1.14.4.3-gompi-2023b

            make ascot5_main CC=h5pcc FLAGS="-foffload=disable"

         For serial version and with IMAS support (under development):

         .. code-block:: bash

            module load zlib/1.3.1-GCCcore-13.3.0 Szip/2.1.1-GCCcore-13.2.0  HDF5/1.14.3-gompi-2023b IMAS-AL-Python/5.4.0-foss-2023b-DD-4.0.0 IPython/8.17.2-GCCcore-13.2.0

            make libascot CC=h5pcc FLAGS="-foffload=disable" MPI=0 -j

            python -m venv --system-site-packages ~/ascot5_venv

            source ~/ascot5_venv/bin/activate

            pip install -e .

            pip install matplotlib pyvista


   .. tab-item:: Lac8 at TCV

      .. card::

         .. code-block:: bash

            make ascot5_main CC=h5cc MPI=0

   .. tab-item:: Marenostrum

      .. card::

         .. code-block:: bash

            module load hdf5/1.8.19 intel/2018.4 impi/2018.4 zlib szip/2.1.1

            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xcommon-avx512 -vecabi=cmdtarget"

   .. tab-item:: Marconi KNL

      .. card::

         .. code-block:: bash

            module load intel/pe-xe-2018--binary intelmpi/2018--binary szip/2.1--gnu--6.1.0 zlib/1.2.8--gnu--6.1.0 hdf5/1.8.18--intelmpi--2018--binary python/3.5.2
            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xmic-avx512 -vecabi=cmdtarget"

   .. tab-item:: Marconi M100 (GPU)

      .. card::

         .. code-block:: bash

            module load xl spectrum_mpi/10.3.1--binary gnu/8.4.0 hdf5/1.12.0--spectrum_mpi--10.3.1--binary szip

   .. tab-item:: Marconi SKL

      .. card::

         With MPI:

         .. code-block:: bash

            module load intel/pe-xe-2020--binary intelmpi/2020--binary gnu/8.3.0 zlib/1.2.11--gnu--8.3.0 szip/2.1.1--gnu--8.3.0 hdf5/1.12.2--intelmpi--2020--binary
            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xcommon-avx512 -vecabi=cmdtarget"

         Without MPI (for working on the login node):

         .. code-block:: bash

            module load load intel/pe-xe-2020--binary gnu/8.3.0 zlib/1.2.11--gnu--8.3.0 szip/2.1.1--gnu--8.3.0 hdf5/1.12.2--intel--pe-xe-2020--binary
            make ascot5_main MPI=0 FLAGS="-qno-openmp-offload -diag-disable 3180"
            make libascot MPI=0 FLAGS="-qno-openmp-offload -diag-disable 3180"

   .. tab-item:: MPCDF Cobra

      .. card::

         .. code-block:: bash

            module load intel/19.1.3 impi/2019.9 git hdf5-mpi
            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180" CC=h5pcc

   .. tab-item:: MPCDF Raven

      .. card::

         .. code-block:: bash

            module load intel/19.1.2 impi/2019.8 git hdf5-mpi anaconda/3/2020.02
            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180"

   .. tab-item:: NERSC Cori

      .. card::

         .. code-block:: bash

            module load cray-hdf5-parallel
            export PMI_NO_FORK=1
            export PMI_NO_PREINITIALIZE=1
            export HDF5_USE_FILE_LOCKING=FALSE
            make ascot5_main CC=h5cc MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180"

   .. tab-item:: OSX (Macports)

      .. card::

         .. code-block:: bash

            port install gcc10
            port install openmpi-gcc10
            port install hdf5 +gcc10 +openmpi +hl

   .. tab-item:: Portal at PPPL

      .. card::

         .. code-block:: bash

            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"
            make libascot MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"

   .. tab-item:: RAT at RFX

      .. card::

         .. code-block:: bash

            module load anaconda
            make ascot5_main MPI=0 CC=h5cc

   .. tab-item:: TOK-cluster at AUG

      .. card::

         .. code-block:: bash

            module load intel/18.0.5 impi/2018.4 hdf5-mpi/1.8.21
            make -j ascot5_main MPI=0 FLAGS="-qno-openmp-offload -diag-disable 3180"

            .. tab-item:: Triton.aalto.fi

         For GCC (outdated):

         .. code-block:: bash

            module load hdf5/1.10.7-openmpi
            make -j ascot5_main MPI=1

         And in the Makefile:

         .. code-block::

            -CFLAGS+=-O2 -lm -Wall -fopenmp -fPIC -std=c11 $(DEFINES) $(FLAGS)
            +CFLAGS+=-O2 -lm -Wall -fopenmp -fPIC -std=c11 -lhdf5_hl $(DEFINES) $(FLAGS)
            -libascot.so: CFLAGS+=-shlib -fPIC -shared
            +libascot.so: CFLAGS+=-fPIC -shared

         For Intel:

         .. code-block:: bash

            module purge
            module load intel-parallel-studio hdf5/1.10.2-openmpi
            export OMPI_MPICC=icc
            make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"

            # Copy binary somewhere safe and use it in batch jobs. These are for the login node
            make clean
            make -j ascot5_main MPI=0 CC=icc FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"
            make -j libascot MPI=0 CC=icc FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"
            module load anaconda

         (Add -xcommon-avx512 to optimize for skl/csl nodes)

   .. tab-item:: EUROfusion gateway

      .. card::

         Serial version:

         For the serail version (without MPI, such as python GUI)

         `module purge ; module load cineca intel/pe-xe-2017--binary intelmpi/2017--binary gnu/6.1.0 zlib/1.2.8--gnu--6.1.0 szip/2.1--gnu--6.1.0 hdf5/1.8.17--gnu--6.1.0 itm-python/3.10`

         .. code-block::

            MPI=0 FLAGS="-I${HDF5_INCLUDE}"

         Parallel version:

         For the parallel version (e.g. to be run on the worker nodes)

         `module purge ; module load cineca intel/pe-xe-2017--binary intelmpi/2017--binary gnu/6.1.0 zlib/1.2.8--gnu--6.1.0 szip/2.1--gnu--6.1.0 hdf5/1.8.17--intelmpi--2017--binary itm-python/3.10`

         .. code-block::

            MPI=1 FLAGS="-I${HDF5_INCLUDE}"

         (this hasn't been really tested, but it is a starting point)

.. _Compilerflags:

Settings when compiling
=======================

Some of the ASCOT5 options require recompiling the code.
Parameters that can be given arguments for ``make`` are (the default values are shown)

.. code-block:: bash

   make -j ascot5_main NSIMD=16 CC=h5cc TARGET=0 VERBOSE=1 MPI=0

.. list-table::
   :widths: 10 50

   * - NSIMD
     - Number of particles simulated in parallel in each SIMD vector.
       These are processed simultaneously by each thread and the optimal number depends on the hardware.
       If unsure, keep the default value.
   * - CC
     - C compiler to use.
   * - TARGET
     - Offload computation to this many Xeon Phi accelerator(s).
       If unsure, do not use this setting.
   * - VERBOSE
     - Print increasing amounts of progress information.

       - 0: No information except bare essentials.
       - 1: Standard information; everything happening outside simulation loops is printed.
       - 2: Extensive information; a record of simulation progress is written to the process-specific \*.stdout file(s).
   * - MPI
     - Enable MPI.
       The code can also be run on multiple nodes without MPI, but doing so requires manual labor.

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
