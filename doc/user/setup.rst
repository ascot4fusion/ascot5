Installing
==========

Installing ASCOT directly from source is recommended for HPC platforms to have all optimizations available.
This is also needed if you want to run the code on GPUs.
Installing via conda is the most convenient option on systems where you only plan to do pre- and postprocessing (e.g. your own laptop) or run only light simulations.

ASCOT has following requirements:

- Python >= 3.12
- Make
- C compiler with OpenMP support (e.g. *gcc*, *icc*, *icx*, *clang*)
- HDF5
- MPI (optional to run multi-node simulations)


From source
***********

Download and compile the code:

.. code-block:: bash

    # <Load any required modules> #
    git clone https://github.com/ascot4fusion/ascot5.git
    cd ascot5
    make ascot

Set up a virtual environment and install the Python package ``a5py`` (which contains ASCOT) along with its dependencies

.. code-block:: bash

    python3 -m venv .ascot-env
    source .ascot-env/bin/activate
    pip install .

For GPUs the process is otherwise the same but you'll need C compiler with CUDA and OpenMP or OpenACC support (e.g. *nvc* for NVIDIA and *clang* for AMD GPUs).
While compiling set ``GPU=1``, and then either enable OpenACC, ``ACC=1``, (for NVIDIA) or use OpenMP, ``OMP=1``, (for AMD).
For example:

.. code-block:: bash

    make ascot GPU=1 ACC=1 CC=nvc

Additional flags such as
   * - ``-xcommon-avx512``, ``-xcore-avx512``, ``-xmic-avx512``
     - Compile the code for Skylake or KNL processors, optimize for Skylake, optimize for KNL.
   * - ``-vecabi=cmdtarget``
     - Enables vector instructions for ``NSIMD`` > 2.
   * - ``-ipo``
     - "Interprocedural Optimization" which might increase the performance somewhat.

are set with ``FLAGS=...``.
For CPUs the number of markers that are vectorized is defined by ``NSIMD`` (default is eight).

Examples of compiling ASCOT can be found (and stored) in the `host-specifics <>`_ repository.

   * - ``-qopt-report=5`` and ``-qopt-report-phase=vec``
     - Generate vectorization reports in \*optrpt files.
       Only useful for developers.


With conda
**********

`Conda <https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html>`_ is a handy package manager that, if not available, can be installed with

If Conda is not available on your cluster, it can be easily `installed <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_ (doesn't require sudo) with:

.. code-block:: bash

    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh

Conda environment for running ASCOT is available on the repository.
To install and activate the conda environment (without MPI):

.. code-block:: bash

    conda env create -f environment.yaml
    conda activate ascot-env

Now you can install ASCOT with ease:

.. code-block:: bash

    make ascot
    pip install .

.. note::
    Add the following lines in your batch script to use the conda environment when submitting jobs via SLURM:

    .. code-block:: bash

       eval "$(conda shell.bash hook)"
       conda activate ascot-env

For MPI it is recommended to use the native MPI implementation instead of one provided by conda, which might work poorly on some platforms.

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


Alternatively one can use MPI implementation packaged with conda.
There is a YAML file readily available in the repository that installs all the necessary packages:

.. code-block:: bash

    conda env create -f environment.yaml
    conda activate ascot-env

To switch to native MPI <https://conda-forge.org/docs/user/tipsandtricks/#using-external-message-passing-interface-mpi-libraries>`_, do the following:



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



Settings when compiling
***********************


Example batch scripts
=====================

.. dropdown:: With MPI

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

.. dropdown:: Without MPI

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


Simulation time
===============

Following table shows some typical run times for ASCOT.

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

Scaling
*******

Due to the fact that i) markers do not interact with each other and ii) the background quantities are not updated based on the marker evolution, ASCOT scales almost perfectly.

.. figure:: ../_static/strong_scaling.png
   :class: with-border
   :align: center
   :width: 500px

   Strong scaling on CPUs and GPUs, see Ref for details.

.. figure:: ../_static/weak_scaling.png
   :class: with-border
   :align: center
   :width: 500px

   Weak scaling on GPUs but same scaling applies to CPUs as well.
   See Ref for details.
   Ref. :cite:t:`tinguely2021modeling`


Parallelization
***************

ASCOT is a high-performance computing code, and understanding its optimization strategies is essential to minimize simulation time and maximize utility.
Given that each marker can be simulated independently, the main optimization strategy is to maximize the number of markers processed concurrently.
The code achieves this through multiple layers of parallelization, as illustrated in Fig. 1.
While the GPU and CPU implementations share similarities, there are key differences in their execution that we will explore, starting with CPU optimizations.

.. grid:: 1 2 2 2

   .. grid-item::

      .. figure:: ../_static/parallelization.svg
         :class: with-border
         :align: center
         :width: 300px

      Fig. 1: How ASCOT5 uses parallelization to simulate multiple markers in parallel.
      The numbers are only for illustration purposes and their values are case specific.

   .. grid-item::

      At the beginning of the simulation, the markers that make up the simulation are initialized.
      The code then spawns *processes* (first parallelization layer) and distributes markers evenly among them.

      Each process places it's markers in a queue and spawns *threads* (second layer).
      A thread fetches a marker from the queue every time it has finished simulating one.

      Each thread is capable of simulate multiple markers concurrently by organizing the markers into a struct of arrays that contain their coordinates.
      This allows the CPU to treat these arrays as *vectors* (third layer), enabling vectorized operations for efficiency.

      The GPU variant starts the same way with *processes* (first layer) being spawned and markers distributed among them.

      However, instead of creating threads immediately and using a queue, all markers assigned to a process are organized into a single large struct of arrays.
      All markers are advanced at the same time, and every time the markers are being operated on, the simulation is *offloaded* to the GPU where *threads* are created (second layer) to perform the operations in parallel.

In the following sections, we will dive deeper into these points and explain their implementation.
For a summary of key points relevant to installation and code usage, jump to the concluding section.

MPI Processes
-------------

In both CPU and GPU variants, the processes are spawned using MPI (Message Passing Interface).
MPI is a protocol whose implementation depends on an external library, commonly ``mpich``, ``openmpi`` or ``intel-mpi``.

ASCOT5 is designed to be a light user of MPI: processes communicate only at the start—before any markers have been simulated—and at the end, when all markers have finished their simulation.
Therefore, the choice of MPI library typically does not impact performance significantly.

What is important to note is that MPI processes are independent and *do not share memory*.
This has implications when using memory-intensive inputs like 3D magnetic fields, 3D wall data, or distribution outputs.
To optimize memory usage, you can run a single MPI process per node, allowing the simulation to access the entire memory available on that node.
However, in some cases, using multiple processes per node may be advantageous, as discussed in the next section.

Another key point is that there is no dynamic load balancing between MPI processes.
For example, if simulating beam ions across three nodes with markers distributed based on beam energy levels (e.g., one node handling 1/3 of the beam energy markers, another handling 1/2, and the last handling full energy markers), the nodes may complete their tasks at different times.
In this scenario, the first node might finish much earlier than the third, but all nodes remain occupied until all tasks are completed.
To mitigate this, it's advisable to randomly permute the markers before writing them to disk or initializing them to ensure a more balanced distribution of work.

Threads (CPU)
-------------

In the CPU variant, the threads are created using OpenMP (Open Multi-Processing).
OpenMP provides directives that serve as "hints" for the compiler on how to parallelize the code, meaning the implementation, or lack of, depends on the compiler.
Compilers commonly used in HPC clusters universally support multithreading.

Unlike MPI processes, threads share memory, and the workload is dynamically balanced through the marker queueing system.
To optimize performance, the number of threads should match the number of physical cores available per processor, assuming a single MPI process is deployed per node (otherwise, set ``Nthreads = Ncores / Nprocesses``).
Some processors support hyperthreading, allowing a single physical core to execute multiple threads in parallel (logical cores).
You can set the number of threads using the ``OMP_NUM_THREADS`` environment variable:

.. code-block:: sh

   export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
   # or to enable hyperthreading
   export OMP_NUM_THREADS=$(( 2 * SLURM_CPUS_PER_TASK ))

.. note::

    In some cases, sharing a node with multiple MPI processes can result in faster execution.
    This could be due to how the cores access the memory or how threads are pinned to cores.
    It is advisable to experiment with different configurations by adjusting the number of processes while keeping the total number of cores (``n x c``) constant.
    This can help determine the optimal setup for your specific simulation.


Vectorization (CPU)
-------------------

During the simulation, a group of markers belonging to the same thread is represented using a struct of arrays, where each array corresponds to a coordinate and the index within the array identifies which marker it belongs to.
The purpose of this struct is to locate markers in contiguous memory, enabling efficient vector operations on them.
These vector operations are executed within SIMD (Single Instruction Multiple Data) loops.

The number of markers within a vector is specified by the compiler flag NSIMD, which should be set to 8 for modern processors that support avx512 instructions.
Note that setting NSIMD=1 effectively disables vectorization.
It is also important to consider that vectorization introduces overhead; thus, one should not expect a linear increase in simulation efficiency with a higher NSIMD value, unlike the behavior observed with threads or processes.

The number of markers within a vector is specified by the compiler flag ``NSIMD``, which should be set to 8 for modern processors that support ``avx512`` instructions.
Note that setting ``NSIMD=1`` effectively disables vectorization.
It is also important to consider that vectorization introduces overhead; thus, one should not expect a linear increase in simulation efficiency with a higher ``NSIMD`` value, unlike the behavior observed with threads or processes.

Threads (GPU)
-------------

In the GPU context, threads function differently from those on the CPU; they are not as independent.
Instead, GPU threads share similarities with vectorized operations, albeit with more flexibility.

Although GPUs allow for the parallel simulation of many more markers than CPUs, computational performance may vary.
Whether the CPU or GPU is faster depends on the available hardware resources and the specifics of the simulation case.
Generally, only large simulations that fully utilize the GPU, i.e. the number of markers is hundreds of thousands at the very least, should be run on it when CPU nodes are also available.

Accelerating simulations
------------------------

1. Choose the appropriate compilation options:

   - MPI: Enable MPI by setting the ``MPI=1`` flag when compiling.
     Ensure that an MPI library such as ``mpich``, ``openmpi``, or ``intel-mpi`` is available.

   - Multithreading: This is enabled by default.
     Make sure to specify the number of threads according to the available cores (see point 2).

   - Vectorization is available if the code is compiled with ``icc`` and the target architechture supports vector operations.
     Set ``NSIMD=8`` and use the compiler flag ``-march=<targetcpu>`` though you may also need to use ``-vecabi=cmdtarget`` or ``-ipo``.

   - GPU (if available): If your cluster supports GPU computing, you can compile for GPU usage by setting ``GPU=1`` and
     either ``OMP=1`` (for OpenMP) or ``ACC=1`` (for OpenACC), depending on what your compiler supports.
     Use GPUs mainly for larger simulations with hundred thousand or more markers.

2. Select the number of processes per node and threads per process for CPU simulations:

   - Ensure that the number of processes multiplied by the number of cores matches the total available cores per node.

   - For memory-intensive simulations, use just one process per node to use the entire node's memory efficiently.

   - Some systems operate more efficiently when the node is shared by multiple processes.
     Test this by varying the number of processes.

   - Test hyperthreading by setting ``OMP_NUM_THREADS`` to twice the number of physical cores.

3. Use number of markers that is a multiple of the number of markers that the system can simulate in parallel.
   Permute the markers randomly before the simulation to ensure balance the work load between processes.


Contributing
============

There is a separate YAML for setting up developer environment:

.. code-block:: bash

    conda env create -f environment-dev.yaml
    conda activate ascot-dev

Within this environment it is possible to run tests and linters, and build the documentation.
The separate developer's guide is built with

.. code-block:: bash

    make doc-dev

and it is located in ``build/doc/dev/index.html`` (it's not available online).
