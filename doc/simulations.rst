.. _Simulations:

===========
Simulations
===========

The simulation is started with:

.. code-block:: bash

   ./ascot5_main --in=in --out=out

or when using ``mpirun`` (or ``srun``):

.. code-block:: bash

   mpirun ./ascot5_main

The default file for reading inputs and storing outputs is ``ascot.h5`` but these can be given explicitly (storing output in the same file with the inputs is recommended):

.. code-block:: bash

   mpirun ./ascot5_main --in=in --out=out

The simulation can be given a description but this can be done in post-processing as well:

.. code-block:: bash

   ascot5_main --d="TEST This is a test run"

The default behavior is to use inputs that are marked as active in the file, but this can be overridden by providing the input QID explicitly:

.. code-block:: bash

   ascot5_main --bfield=0123456789

When not using MPI, it is possible to divide the simulation into separate processes manually with:

.. code-block:: bash

   ascot5_main --mpi_size=size --mpi_rank=rank

.. _Simulationoptions:

Options
=======

Simulation options are listed below.
Click on the parameter name for more details.
Note that the first underscore does not belong in a parameter name and only appear here for technical reasons.

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

.. _Batchjobs:

Batch jobs
==========

Use this script to run ASCOT5 in multiple nodes using MPI (MPI=1 when the code was compiled).

.. code-block::

   #!/bin/bash

   ## Set required nodes and CPUs (number of processes == number of nodes)
   #SBATCH -N10 -n10 -c48

   #SBATCH -t 24:00:00
   #SBATCH -J ascot5

   #SBATCH -e %x.e%j
   #SBATCH -o %x.o%j

   export <your exports>
   module load <your modules>

   echo Job name $SLURM_JOB_NAME
   echo Job id $SLURM_JOB_ID

   # Some platforms require that the number of threads must be given explicitly
   # In those cases use (1-2) x number of CPUs (check which is faster)
   export OMP_NUM_THREADS=48

   INPUTFILE=ascot

   date
   export FOR_PRINT=$SLURM_JOB_ID.stdout
   mpirun ./ascot5_main --in=$INPUTFILE --d="YOURTAG Your description"
   date

Use this script to run ASCOT5 in multiple nodes without MPI (MPI=0 when the code was compiled).
This results in multiple output files that you have to combine using the python tool ``a5combine``.

.. code-block::

   #!/bin/bash

   ## How many jobs this run is divided into
   NJOBS=100

   ## Name of the input file. Each output file begins with this name followed by
   ## its "mpi_rank"
   INPUTFILE=ascot

   for i in $(seq 0 $(($NJOBS-1)))
   do
   filename=${RANDOM}
   cat > $filename << EOF
   #!/bin/bash

   ## Use only single node and single process
   #SBATCH -N1 -n1 -c48

   #SBATCH -t 24:00:00
   #SBATCH -J ascot5

   #SBATCH -e ${i}.e%j
   #SBATCH -o ${i}.o%j

   export <your exports>
   module load <your modules>

   # Some platforms require that the number of threads must be given explicitly
   # In those cases use (1-2) x number of CPUs (check which is faster)
   export OMP_NUM_THREADS=96

   echo Job name $SLURM_JOB_NAME
   echo Job id $SLURM_JOB_ID

   date
   ./ascot5_main --mpi_rank=$i --mpi_size=$NJOBS --in=$INPUTFILE --d="YOURTAG Your description"
   date

   EOF

   sbatch $filename
   rm $filename
   echo $i
   ## If multiple jobs are reading from the same file at the same time this can
   ## cause problems in the system, so give system some time to breathe
   sleep 5s
   done

.. _ExampleRuntimes:

Examples
========

.. list-table:: Examples of simulation times and number of markers used
   :widths: 50 25 25 25
   :header-rows: 1

   * - Simulation
     - Hardware
     - Number of markers
     - CPU time
   * - Alpha particle slowing-down in 3D ITER (Not real numbers!)
     - Xeon-Phi, 10 nodes (MARCONI)
     - :math:`1\times10^6`
     - 24 h
   * -
     -
     -
     -
