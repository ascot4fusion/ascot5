.. _Installing:

==========
Installing
==========

1. Clone the repository.

   .. code-block:: bash

      git clone git@version.aalto.fi:ascot/ascot5.git ascot5
      cd ascot5

2. Compile the main simulation program ``ascot5_main`` and the library ``libascot.so`` (which is needed for post-processing so it is not needed to execute simulations).

   .. code-block:: bash

      make ascot5_main
      make libascot

   Compiling requires **C compiler**, **HDF5**, **OpenMP** and **MPI** (optional), and for the postprocessing one needs **Python 3** and **VTK** (optional).
   See :ref:`here<Compiling>` for tips on how to compile the code (on different platforms).

3. Everytime there is a new release, you can update the code as:

   .. code-block:: bash

      git pull
      git checkout master
      make clean
      make ascot5_main
      make libascot

   Always use the branch ``master`` when running simulations unless you specifically need something from a feature branch.

4. Set up ``a5py`` Python package needed for input generation and post-processing.

   .. code-block:: bash

      cd ..

   It is recommended to set a virtual environment when working with ASCOT5. You can use either ``virtualenv`` (3rd party tool)

   .. code-block:: bash

      virtualenv -p python3 --system-site-packages ascotenv
      source ascotenv/bin/activate

   or ``venv`` (packs with Python 3)

   .. code-block:: bash

      python3 -m venv --system-site-packages ascotenv
      source ascotenv/bin/activate

   or Conda (preferred)

   Edit ``ascot5/a5py/ascot5.yaml`` file

   .. code-block:: bash

      conda create --file ascot5/a5py/requirements.txt --name ascotenv
      conda activate ascotenv

   Whichever was chosen, the python package can be installed with

   .. code-block:: bash

      pip install -e ascot5/a5py

5. Set environment variables.

   Note the directory where ``libascot.so`` is located as it needs to be added to ``LD_LIBRARY_PATH``:

   .. code-block:: bash

      cd ascot5
      echo $PWD

   Edit your `.bashrc` and/or `.bash_profile`

   .. code-block:: bash

      emacs -nw ~/.bashrc

   Add the following lines (along with any ``module load`` or ``export`` that was needed for the code to compile

   .. code-block::

      <module loads here>
      export LD_LIBRARY_PATH=/path/to/ascot5folder:$LD_LIBRARY_PATH
      export EDITOR=/usr/bin/emacs

   Last line is optional and it sets Emacs as the editor when editing ASCOT5 options. You can also add

   .. code-block::

      source activate /path/to/ascot5env # If using venv or virtualenv
      conda activate ascot5env # If using Conda

   which automatically activates ASCOT5 environment each time you open the terminal or login (otherwise you have to execute the line manually).

   Close the terminal (or log out) and the environment is set when you open a new one (or log in).

6. Test that ASCOT5 is working by running the tutorial.

.. _Compiling:

================================
Compiling on different platforms
================================

Add your platform here!

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

=======================
Settings when compiling
=======================

Some of the ASCOT5 options require recompiling the code.
Parameters that can be given arguments for ``make`` are

.. code-block:: bash

   make -j ascot5_main NSIMD=16 CC=icc TARGET=1 VERBOSE=1 MPI=1 NOGIT=1

.. list-table::
   :widths: 10 50

   * - NSIMD
     - Number of particles in a group. These are processed simultaneously by each thread and the optimal number depends on the platform. If unsure, keep the default value.
   * - CC
     - Compiler.
   * - TARGET
     - Offload computation to Xeon Phi accelerator.
   * - VERBOSE
     - Print increasing amounts of progress information. 0: No information except bare essentials. 1: Standard information; everything happening outside simulation loops is printed. 2: Extensive information; a record of simulation progress is written process-specific \*.stdout file(s).
   * - MPI
     - Enable MPI.
   * - NOGIT
     - Disable recording of repository status if Git is not available.

Additional parameters can be found in ``ascot5.h``, but there is rarely a need to change these.

.. doxygendefine:: MAX_SPECIES

.. doxygendefine:: MHD_MODES_MAX_NUM

.. doxygendefine:: WIENERSLOTS

.. doxygendefine:: A5_CCOL_USE_GEOBM

.. doxygendefine:: A5_EXTREMELY_SMALL_TIMESTEP

.. doxygendefine:: A5_PRINTPROGRESSINTERVAL

.. doxygendefine:: A5_WTIME

.. doxygendefine:: INTERP_SPL_EXPL

.. doxygendefine:: A5_CCOL_USE_TABULATED

.. _Simulationoptions:

=======
Options
=======

.. autoproperty:: a5py.ascot5io.options.Opt.SIM_MODE
   :noindex:

.. _Simulations:

===========
Simulations
===========

.. _Batchjobs:

==========
Batch jobs
==========
