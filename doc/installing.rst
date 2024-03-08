.. _Installing:

Installation
============

Full installation
*****************

The most convenient way to install ASCOT5 is to use Conda:

.. code-block:: bash

   git clone https://github.com/ascot4fusion/ascot5.git
   cd ascot5
   conda env create -f environment.yaml
   conda activate ascot-env
   make libascot -j
   make ascot5_main -j # Also bbnbi5 if needed
   pip install -e .

To use the Conda environment in your batch script, add the following lines:

.. code-block:: bash

   eval "$(conda shell.bash hook)"
   conda activate ascot-env

The simulation options are edited with the local text editor, which usually happens to be ``vim``.
Consider adding the following line to your ``.bashrc`` (or ``.bash_profile`` if working locally):

.. code-block:: bash

   export EDITOR=/usr/bin/emacs

.. note::
   Whenever there is a new release, you can update the code as:

   .. code-block:: bash

      git pull
      make clean
      make ascot5_main -j
      make libascot -j

   Always use the ``main`` branch when running simulations unless you specifically need something from a feature branch.

Test that ASCOT5 was properly installed by running the :ref:`introduction<Tutorial>`.

.. note::
   The default instructions install ASCOT5 without MPI, but MPI is recommended for running simulations on multiple nodes.

   You can use either MPI packaged with Conda (``openmpi`` or ``mpich``):

   .. code-block:: bash

      conda install hdf5=*=*openmpi* mpi=*=*openmpi*

   or the native MPI library:

   .. code-block:: bash

      module load <local-mpi-library>
      mpirun -V # Prints the version of the MPI library to be used in the next line
      conda install "<local-mpi-library>=x.y.z=external_*"
      which mpicc # Copy this path to the following commands
      export HDF5_CC=path_to_mpicc
      export HDF5_CLINKER=path_to_mpicc

   Conda does not have all versions of MPI libraries available, so one might have to use asterix in minor version number, e.g.

   .. code-block:: bash

      conda install "openmpi=4.*=external_*"

.. admonition:: Troubleshooting

   Conda complains that ``freeqdsk`` was not found or could not be installed.

   - Remove ``freeqdsk`` from ``environment.yaml`` and install it with ``pip`` once the environment is activated.
     The cause of this error is unknown.

   Make stops since ``libhdf5*.a`` could not be located.

   - The compiler should be using shared libraries, not static.
     Add the following flag for the compiler: ``make ascot5_main FLAGS="-shlib"``

Minimal installation (for HPC usage)
************************************

A typical way to use ASCOT5 is to run simulations on a computing cluster and carry out pre- and postprocessing on a home cluster or a workstation.
Building the code from the source (without Conda) provides a light installation with the capability to run simulations using native libraries and perform limited pre- and post-processing.

.. rubric:: Requirements

- C compiler
- HDF5
- OpenMP
- MPI
- Python3.10

1. Install the requirements or use the module system.

2. Download the source and set up the virtual environment:

   .. code-block:: bash

      git clone https://github.com/ascot4fusion/ascot5.git
      python -m venv ascot-env
      source activate ascot-env/bin/activate

3. Install ``a5py`` and compile the executables which will be located at ``build/``:

   .. code-block:: bash

      cd ascot5
      pip install -e .
      make ascot5_main -j # Also bbnbi5 if needed

See :ref:`here<Compiling>` for tips on how to compile the code on different platforms.

(**Optional**) Add the following lines to your `.bashrc`:

.. code-block::

   <module loads and exports here>
   source activate /path/to/ascot5env

This will automatically activate ASCOT5 environment each time you open a terminal or login.


For Developers
**************

This is a full installation from the source (using Conda) and with the optional packages present that are required to build ``ascot2py.py`` and the documentation.
Note that the first step requires you to add SSH keys on GitHub whenever on a new machine.

.. code-block:: bash

   git clone git@github.com:ascot4fusion/ascot5.git
   cd ascot5
   conda env create -f environment-dev.yaml
   conda activate ascot-dev
   make libascot -j
   make ascot5_main -j
   pip install -e .

.. _Compiling:

Compiling on different platforms
================================

ASCOT5 doesn't support CMake (yet) so it is up for the user to provide the required libraries.
Here we have listed some platforms where ASCOT5 has been used and how it was compiled.

CSC.fi puhti
************

.. code-block:: bash

   module load StdEnv intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4 hdf5/1.10.4-mpi python-data

   make -j ascot5_main MPI=1

Alternatively:

.. code-block:: bash

   make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"

CSD3
****

.. code-block:: bash

   module load hdf5/mpi/intel/2019.3/1.10.5
   make ascot5_main MPI=1 CC=mpicc FLAGS="-no-multibyte-chars -qno-openmp-offload -diag-disable 3180 -xmic-avx512 -vecabi=cmdtarget" LFLAGS="-lhdf5_hl -lhdf5" -j
   make libascot MPI=1 CC=mpicc FLAGS="-no-multibyte-chars -qno-openmp-offload -diag-disable 3180" LFLAGS="-lhdf5_hl -lhdf5" -j

Freia (UKAEA) (work in progress)
********************************

.. code-block:: bash

   make ascot5_main MPI=0 CC=gcc

For libascot, one user needed to revert to python/3.5.1 and command

.. code-block:: bash

   make libascot MPI=0

ITER sdcc
*********

Going with intel compilers.

.. code-block:: bash

   module load GCCcore/11.3.0 zlib/1.2.12-GCCcore-11.3.0 binutils/2.38-GCCcore-11.3.0 intel-compilers/2022.1.0 numactl/2.0.14-GCCcore-11.3.0 UCX/1.12.1-GCCcore-11.3.0 impi/2021.6.0-intel-compilers-2022.1.0 iimpi/2022a Szip/2.1.1-GCCcore-11.3.0 HDF5/1.13.1-iimpi-2022a

   make ascot5_main CC=h5pcc FLAGS=-qno-openmp-offload -diag-disable 3180

Lac8 at TCV
***********

.. code-block:: bash

   make ascot5_main CC=h5cc MPI=0

Marenostrum (WIP)
*****************

.. code-block:: bash

   module load hdf5/1.8.19 intel/2018.4 impi/2018.4 zlib szip/2.1.1

   make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xcommon-avx512 -vecabi=cmdtarget"

Marconi KNL
***********

.. code-block:: bash

   module load intel/pe-xe-2018--binary intelmpi/2018--binary szip/2.1--gnu--6.1.0 zlib/1.2.8--gnu--6.1.0 hdf5/1.8.18--intelmpi--2018--binary python/3.5.2
   make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xmic-avx512 -vecabi=cmdtarget"

Marconi M100 (GPU)
******************

.. code-block:: bash

   module load xl spectrum_mpi/10.3.1--binary gnu/8.4.0 hdf5/1.12.0--spectrum_mpi--10.3.1--binary szip

Marconi SKL
***********

With MPI:

.. code-block:: bash

   module load intel/pe-xe-2020--binary intelmpi/2020--binary gnu/8.3.0 zlib/1.2.11--gnu--8.3.0 szip/2.1.1--gnu--8.3.0 hdf5/1.12.2--intelmpi--2020--binary
   make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xcommon-avx512 -vecabi=cmdtarget"

Without MPI (for working on the login node):

.. code-block:: bash

   module load load intel/pe-xe-2020--binary gnu/8.3.0 zlib/1.2.11--gnu--8.3.0 szip/2.1.1--gnu--8.3.0 hdf5/1.12.2--intel--pe-xe-2020--binary
   make ascot5_main MPI=0 FLAGS="-qno-openmp-offload -diag-disable 3180"
   make libascot MPI=0 FLAGS="-qno-openmp-offload -diag-disable 3180"

MPCDF Cobra
***********

.. code-block:: bash

   module load intel/19.1.3 impi/2019.9 git hdf5-mpi
   make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180" CC=h5pcc

MPCDF Raven
***********

.. code-block:: bash

   module load intel/19.1.2 impi/2019.8 git hdf5-mpi anaconda/3/2020.02
   make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180"

NERSC Cori
**********

.. code-block:: bash

   module load cray-hdf5-parallel
   export PMI_NO_FORK=1
   export PMI_NO_PREINITIALIZE=1
   export HDF5_USE_FILE_LOCKING=FALSE
   make ascot5_main CC=h5cc MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180"

OSX (Macports)
**************

.. code-block:: bash

   port install gcc10
   port install openmpi-gcc10
   port install hdf5 +gcc10 +openmpi +hl

Portal at PPPL
**************

.. code-block:: bash

   make ascot5_main MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"
   make libascot MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"

RAT at RFX
**********

.. code-block:: bash

   module load anaconda
   make ascot5_main MPI=0 CC=h5cc

TOK-cluster at AUG
******************

.. code-block:: bash

   module load intel/18.0.5 impi/2018.4 hdf5-mpi/1.8.21
   make -j ascot5_main MPI=0 FLAGS="-qno-openmp-offload -diag-disable 3180"

Triton.aalto.fi
***************

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

.. _Compilerflags:

EUROfusion gateway
******************

Serial version
--------------
For the serail version (without MPI, such as python GUI)

`module purge ; module load cineca intel/pe-xe-2017--binary intelmpi/2017--binary gnu/6.1.0 zlib/1.2.8--gnu--6.1.0 szip/2.1--gnu--6.1.0 hdf5/1.8.17--gnu--6.1.0 itm-python/3.10`

.. code-block::

   MPI=0 FLAGS="-I${HDF5_INCLUDE}"

Parallel version
----------------
For the parallel version (e.g. to be run on the worker nodes)

`module purge ; module load cineca intel/pe-xe-2017--binary intelmpi/2017--binary gnu/6.1.0 zlib/1.2.8--gnu--6.1.0 szip/2.1--gnu--6.1.0 hdf5/1.8.17--intelmpi--2017--binary itm-python/3.10`

.. code-block::

   MPI=1 FLAGS="-I${HDF5_INCLUDE}"

(this hasn't been really tested, but it is a starting point)

Settings when compiling
=======================

Some of the ASCOT5 options require recompiling the code.
Parameters that can be given arguments for ``make`` are (the default values are shown)

.. code-block:: bash

   make -j ascot5_main NSIMD=16 CC=icc TARGET=0 VERBOSE=1 MPI=1 NOGIT=0

.. list-table::
   :widths: 10 50

   * - NSIMD
     - Number of particles simulated in parallel in each SIMD vector.
       These are processed simultaneously by each thread and the optimal number depends on the hardware.
       If unsure, keep the default value.
   * - CC
     - The compiler.
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
       The code can be run on multiple nodes without MPI, but doing so requires manual labor.
   * - NOGIT
     - Disable recording of repository status if Git is not available.

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

Additional compile-time parameters can be found in ``ascot5.h``, but there is rarely a need to change these.

.. doxygendefine:: MAX_SPECIES

.. doxygendefine:: MHD_MODES_MAX_NUM

.. doxygendefine:: WIENERSLOTS

.. doxygendefine:: A5_CCOL_USE_GEOBM

.. doxygendefine:: A5_EXTREMELY_SMALL_TIMESTEP

.. doxygendefine:: A5_PRINTPROGRESSINTERVAL

.. doxygendefine:: A5_WTIME

.. doxygendefine:: INTERP_SPL_EXPL

.. doxygendefine:: A5_CCOL_USE_TABULATED
