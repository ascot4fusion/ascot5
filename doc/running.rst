.. _Installing:

==========
Installing
==========

1. Clone the repository.

   .. code-block:: bash

      git clone git@github.com:ascot4fusion/ascot5.git ascot5
      cd ascot5

2. Compile the main simulation program ``ascot5_main`` and the library ``libascot.so`` (which is needed for post-processing but not for executing simulations).

   .. code-block:: bash

      make ascot5_main
      make libascot

   Compiling requires **C compiler**, **HDF5**, **OpenMP** and **MPI** (optional).
   For postprocessing one needs **Python 3** and **VTK** (optional for 3D walls).
   See :ref:`here<Compiling>` for tips on how to compile the code on different platforms.
   The binaries will be located in ``build`` folder.

3. Whenever there is a new release, you can update the code as:

   .. code-block:: bash

      git pull
      git checkout main
      make clean
      make ascot5_main
      make libascot

   Always use the branch ``main`` when running simulations unless you specifically need something from a feature branch.

4. Set up a virtual environment

   .. code-block:: bash

      cd ..

   You can use either ``virtualenv`` (3rd party tool)

   .. code-block:: bash

      virtualenv -p python3 --system-site-packages ascotenv
      source ascotenv/bin/activate

   or ``venv`` (packs with Python 3)

   .. code-block:: bash

      python3 -m venv --system-site-packages ascotenv
      source ascotenv/bin/activate

   or Conda

   .. code-block:: bash

      conda create --name ascotenv
      conda activate ascotenv

5. Install ``a5py`` Python package needed for input generation and post-processing.

   .. code-block:: bash

      pip install -e ascot5

   For developers, the extra packages required for generating the documentation are installed with

   .. code-block:: bash

      pip install -e ascot5[doc]

   In addition Doxygen is needed.

6. (Optional) Edit your `.bashrc` and/or `.bash_profile`

   .. code-block:: bash

      emacs -nw ~/.bashrc

   Add the following lines (along with any ``module load`` or ``export`` that was needed for the code to compile

   .. code-block::

      <module loads and exports here>
      export EDITOR=/usr/bin/emacs # To use emacs when editing options
      source activate /path/to/ascot5env # If using venv or virtualenv
      conda activate ascot5env # If using Conda

   This will automatically activate ASCOT5 environment each time you open a terminal or login.

   Close the terminal (or log out) and the environment is set when you open a new one (or log in).

7. Test that ASCOT5 is working by running the :ref:`tutorial<Tutorial>`.

.. _Compiling:

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

==================
Simulation options
==================

.. currentmodule:: a5py.ascot5io.options

.. rubric:: Simulation mode and time step

.. autosummary::

   ~Opt._SIM_MODE
   ~Opt._ENABLE_ADAPTIVE
   ~Opt._RECORD_MODE
   ~Opt._FIXEDSTEP_USE_USERDEFINED
   ~Opt._FIXEDSTEP_USERDEFINED
   ~Opt._FIXEDSTEP_GYRODEFINED
   ~Opt._ADAPTIVE_TOL_ORBIT
   ~Opt._ADAPTIVE_TOL_CCOL
   ~Opt._ADAPTIVE_MAX_DRHO
   ~Opt._ADAPTIVE_MAX_DPHI

.. rubric:: Simulation end conditions

.. autosummary::

   ~Opt._ENDCOND_SIMTIMELIM
   ~Opt._ENDCOND_CPUTIMELIM
   ~Opt._ENDCOND_RHOLIM
   ~Opt._ENDCOND_ENERGYLIM
   ~Opt._ENDCOND_WALLHIT
   ~Opt._ENDCOND_MAXORBS
   ~Opt._ENDCOND_NEUTRALIZED
   ~Opt._ENDCOND_IONIZED
   ~Opt._ENDCOND_LIM_SIMTIME
   ~Opt._ENDCOND_MAX_MILEAGE
   ~Opt._ENDCOND_MAX_CPUTIME
   ~Opt._ENDCOND_MAX_RHO
   ~Opt._ENDCOND_MIN_RHO
   ~Opt._ENDCOND_MIN_ENERGY
   ~Opt._ENDCOND_MIN_THERMAL
   ~Opt._ENDCOND_MAX_TOROIDALORBS
   ~Opt._ENDCOND_MAX_POLOIDALORBS

.. rubric:: Active physics

.. autosummary::

   ~Opt._ENABLE_ORBIT_FOLLOWING
   ~Opt._ENABLE_COULOMB_COLLISIONS
   ~Opt._ENABLE_MHD
   ~Opt._ENABLE_ATOMIC
   ~Opt._DISABLE_FIRSTORDER_GCTRANS
   ~Opt._DISABLE_ENERGY_CCOLL
   ~Opt._DISABLE_PITCH_CCOLL
   ~Opt._DISABLE_GCDIFF_CCOLL

.. rubric:: Distributions

.. autosummary::

   ~Opt._ENABLE_DIST_5D
   ~Opt._ENABLE_DIST_6D
   ~Opt._ENABLE_DIST_RHO5D
   ~Opt._ENABLE_DIST_RHO6D
   ~Opt._DIST_MIN_R
   ~Opt._DIST_MAX_R
   ~Opt._DIST_NBIN_R
   ~Opt._DIST_MIN_PHI
   ~Opt._DIST_MAX_PHI
   ~Opt._DIST_NBIN_PHI
   ~Opt._DIST_MIN_Z
   ~Opt._DIST_MAX_Z
   ~Opt._DIST_NBIN_Z
   ~Opt._DIST_MIN_RHO
   ~Opt._DIST_MAX_RHO
   ~Opt._DIST_NBIN_RHO
   ~Opt._DIST_MIN_THETA
   ~Opt._DIST_MAX_THETA
   ~Opt._DIST_NBIN_THETA
   ~Opt._DIST_MIN_PPA
   ~Opt._DIST_MAX_PPA
   ~Opt._DIST_NBIN_PPA
   ~Opt._DIST_MIN_PPE
   ~Opt._DIST_MAX_PPE
   ~Opt._DIST_NBIN_PPE
   ~Opt._DIST_MIN_PR
   ~Opt._DIST_MAX_PR
   ~Opt._DIST_NBIN_PR
   ~Opt._DIST_MIN_PPHI
   ~Opt._DIST_MAX_PPHI
   ~Opt._DIST_NBIN_PPHI
   ~Opt._DIST_MIN_PZ
   ~Opt._DIST_MAX_PZ
   ~Opt._DIST_NBIN_PZ
   ~Opt._DIST_MIN_TIME
   ~Opt._DIST_MAX_TIME
   ~Opt._DIST_NBIN_TIME
   ~Opt._DIST_MIN_CHARGE
   ~Opt._DIST_MAX_CHARGE
   ~Opt._DIST_NBIN_CHARGE

.. rubric:: Constant of motion distribution

.. autosummary::

   ~Opt._ENABLE_DIST_COM
   ~Opt._DIST_MIN_MU
   ~Opt._DIST_MAX_MU
   ~Opt._DIST_NBIN_MU
   ~Opt._DIST_MIN_EKIN
   ~Opt._DIST_MAX_EKIN
   ~Opt._DIST_NBIN_EKIN
   ~Opt._DIST_MIN_PTOR
   ~Opt._DIST_MAX_PTOR
   ~Opt._DIST_NBIN_PTOR

.. rubric:: Recording marker trajectories

.. autosummary::

   ~Opt._ENABLE_ORBITWRITE
   ~Opt._ORBITWRITE_MODE
   ~Opt._ORBITWRITE_NPOINT
   ~Opt._ORBITWRITE_POLOIDALANGLES
   ~Opt._ORBITWRITE_TOROIDALANGLES
   ~Opt._ORBITWRITE_RADIALDISTANCES
   ~Opt._ORBITWRITE_INTERVAL

.. rubric:: Transport coefficients

.. autosummary::

   ~Opt._ENABLE_TRANSCOEF
   ~Opt._TRANSCOEF_INTERVAL
   ~Opt._TRANSCOEF_NAVG
   ~Opt._TRANSCOEF_RECORDRHO

.. currentmodule:: a5py

.. _Simulations:

===========
Simulations
===========

TBD

.. _Batchjobs:

==========
Batch jobs
==========

TBD
