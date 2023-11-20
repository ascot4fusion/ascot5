.. _Installing:

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

See :ref:`here<Compiling>` for tips on how to compile the code on different platforms.

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

.. note::
   Whenever there is a new release, you can update the code as:

   .. code-block:: bash

      cd ascot5
      git pull
      git checkout main
      make clean
      make ascot5_main
      make libascot

   Always use the branch ``main`` when running simulations unless you specifically need something from a feature branch.

.. rubric:: (Optional)

Add the following lines to your `.bashrc` and/or `.bash_profile`:

.. code-block::

   <module loads and exports here>
   export EDITOR=/usr/bin/emacs       # To use emacs when editing options
   source activate /path/to/ascot5env

This will automatically activate ASCOT5 environment each time you open a terminal or login.

Test that ASCOT5 is working by running the :ref:`introduction<Tutorial>`.

.. _Compiling:

Compiling on different platforms
================================

ASCOT5 doesn't support CMake (yet) so it is up for the user to provide the required libraries.
Here we have listed some platforms where ASCOT5 has been used and how it was compiled.

Aalto desktops
**************

.. code-block:: bash

   pkcon install hdf5-helpers
   pkcon install libhdf5-dev
   pkcon install hdf5-tools
   make -j ascot5_main VERBOSE=1 FLAGS="-foffload=disable"

CSC.fi puhti
************

.. code-block:: bash

   module load StdEnv intel/19.0.4  hpcx-mpi/2.4.0  intel-mkl/2019.0.4  hdf5/1.10.4-mpi python-data

   make -j ascot5_main MPI=1 VERBOSE=1

Alternatively:

.. code-block:: bash

   make ascot5_main VERBOSE=2 MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"

Freia (UKAEA) (work in progress)
********************************

.. code-block:: bash

   make ascot5_main MPI=0 VERBOSE=1 CC=gcc

For libascot, one user needed to revert to python/3.5.1 and command

.. code-block:: bash

   make libascot MPI=0 VERBOSE=1

ITER sdcc
*********

Going with intel compilers.

.. code-block:: bash

   module load GCCcore/11.3.0 zlib/1.2.12-GCCcore-11.3.0 binutils/2.38-GCCcore-11.3.0 intel-compilers/2022.1.0 numactl/2.0.14-GCCcore-11.3.0 UCX/1.12.1-GCCcore-11.3.0 impi/2021.6.0-intel-compilers-2022.1.0 iimpi/2022a Szip/2.1.1-GCCcore-11.3.0 HDF5/1.13.1-iimpi-2022a

   make ascot5_main CC=h5pcc FLAGS=-qno-openmp-offload -diag-disable 3180

Lac8 at TCV
***********

.. code-block:: bash

   make ascot5_main CC=h5cc MPI=0 VERBOSE=2

Marenostrum (WIP)
*****************

.. code-block:: bash

   module load hdf5/1.8.19 intel/2018.4 impi/2018.4 zlib szip/2.1.1

   make ascot5_main VERBOSE=2 MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xcommon-avx512 -vecabi=cmdtarget"

Marconi KNL
***********

.. code-block:: bash

   module load intel/pe-xe-2018--binary intelmpi/2018--binary szip/2.1--gnu--6.1.0 zlib/1.2.8--gnu--6.1.0 hdf5/1.8.18--intelmpi--2018--binary python/3.5.2
   make ascot5_main VERBOSE=2 MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xmic-avx512 -vecabi=cmdtarget"

Marconi M100 (GPU)
******************

.. code-block:: bash

   module load xl spectrum_mpi/10.3.1--binary gnu/8.4.0 hdf5/1.12.0--spectrum_mpi--10.3.1--binary szip

Marconi SKL
***********

.. code-block:: bash

   module load intel/pe-xe-2020--binary intelmpi/2020--binary gnu/8.3.0 zlib/1.2.11--gnu--8.3.0 szip/2.1.1--gnu--8.3.0  hdf5/1.12.2--intelmpi--2020--binary python/3.5.2
   make ascot5_main VERBOSE=2 MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -xcommon-avx512 -vecabi=cmdtarget"

MPCDF Cobra
***********

.. code-block:: bash

   module load intel/19.1.3 impi/2019.9 git hdf5-mpi
   make ascot5_main MPI=1 VERBOSE=1 FLAGS="-qno-openmp-offload -diag-disable 3180" CC=h5pcc

MPCDF Raven
***********

.. code-block:: bash

   module load intel/19.1.2 impi/2019.8 git hdf5-mpi anaconda/3/2020.02
   make ascot5_main MPI=1 VERBOSE=1 FLAGS="-qno-openmp-offload -diag-disable 3180"

NERSC Cori
**********

.. code-block:: bash

   module load cray-hdf5-parallel
   export PMI_NO_FORK=1
   export PMI_NO_PREINITIALIZE=1
   export HDF5_USE_FILE_LOCKING=FALSE
   make ascot5_main CC=h5cc VERBOSE=0 MPI=1 FLAGS="-qno-openmp-offload –diag-disable 3180"

OSX (Macports)
**************

.. code-block:: bash

   port install gcc10
   port install openmpi-gcc10
   port install hdf5 +gcc10 +openmpi +hl

RAT at RFX
**********

.. code-block:: bash

   module load anaconda
   make ascot5_main MPI=0 VERBOSE=1 CC=h5cc

TOK-cluster at AUG
******************

.. code-block:: bash

   module load intel/18.0.5 impi/2018.4 hdf5-mpi/1.8.21
   make -j ascot5_main MPI=0 VERBOSE=1 FLAGS="-qno-openmp-offload -diag-disable 3180"

Triton.aalto.fi
***************

For GCC:

.. code-block:: bash

   module load hdf5/1.10.7-openmpi
   make -j ascot5_main MPI=1 VERBOSE=1

Makefile:

.. code-block::

   -CFLAGS+=-O2 -lm -Wall -fopenmp -fPIC -std=c11 $(DEFINES) $(FLAGS)
   +CFLAGS+=-O2 -lm -Wall -fopenmp -fPIC -std=c11 -lhdf5_hl $(DEFINES) $(FLAGS)
   -libascot.so: CFLAGS+=-shlib -fPIC -shared
   +libascot.so: CFLAGS+=-fPIC -shared

For Intel:

.. code-block:: bash

   module load intel-parallel-studio hdf5/1.10.2-openmpi
   export OMPI_MPICC=icc
   make ascot5_main VERBOSE=2 MPI=1 FLAGS="-qno-openmp-offload -diag-disable 3180 -vecabi=cmdtarget"

(add -xcommon-avx512 to optimize for skl/csl nodes)

Vdiubuntu.aalto.fi
******************

Compiling libascot.so requires that you change Makefile as

.. code-block::

   libascot.so: libascot.o $(OBJS)
   - $(CC) $(CFLAGS) -o $@ $^
   + $(CC) $(CFLAGS) -o $@ $^ -lhdf5_hl -lhdf5


.. code-block:: bash

   module load hdf5
   make VERBOSE=2 MPI=0 libascot.so CC=gcc FLAGS="-foffload=disable" -j

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

Compiler flags can be provided with ``FLAGS`` parameter, e.g.

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
