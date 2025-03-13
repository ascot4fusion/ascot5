============
Optimization
============

ASCOT5 is a high-performance computing (HPC) code, and understanding its optimization strategies is essential to minimize simulation time and maximize efficiency.
Given that each marker can be simulated independently, the main optimization strategy is to maximize the number of markers processed concurrently.
The code achieves this through multiple layers of parallelization, as illustrated in Fig. 1.
While the GPU and CPU implementations share similarities, there are key differences in their execution that we will explore, starting with CPU optimizations.

.. grid:: 1 2 2 2

   .. grid-item::

      .. figure:: ../figures/ascot5parallel.svg
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
*************

In both CPU and GPU variants, the processes are spawned using MPI (Message Passing Interface).
MPI is a protocol whose implementation depends on an external library, commonly ``mpich``, ``openmpi`` or ``intel-mpi``.
To enable MPI, ensure that an MPI library is available and set the ``MPI=1`` flag during compilation.

ASCOT5 is designed to be a light user of MPI: processes communicate only at the start—before any markers have been simulated—and at the end, when all markers have finished their simulation.
Therefore, the choice of MPI library typically does not impact performance significantly.

.. admonition:: Glossary

   **mpirun/mpiexec:** A program used to launch MPI processes.

   **mpicc:** A wrapper compiler that links the source code with the MPI library.

What is important to note is that MPI processes are independent and *do not share memory*.
This has implications when using memory-intensive inputs like 3D magnetic fields, 3D wall data, or distribution outputs.
To optimize memory usage, you can run a single MPI process per node, allowing the simulation to access the entire memory available on that node.
However, in some cases, using multiple processes per node may be advantageous, as discussed in the next section.

Another key point is that there is no dynamic load balancing between MPI processes.
For example, if simulating beam ions across three nodes with markers distributed based on beam energy levels (e.g., one node handling 1/3 of the beam energy markers, another handling 1/2, and the last handling full energy markers), the nodes may complete their tasks at different times.
In this scenario, the first node might finish much earlier than the third, but all nodes remain occupied until all tasks are completed.
To mitigate this, it's advisable to randomly permute the markers before writing them to disk or initializing them to ensure a more balanced distribution of work.

*For Developers:* When using ``ascot5_main``, each MPI process reads its input files separately, including marker data.
The MPI rank only determines which markers are assigned for the simulation.
The C kernel uses MPI to collect all simulation results, such as marker end states, and gather them to the root process.

To launch simulations using MPI via the Python interface, users need a basic understanding of MPI commands:

.. tab-set::

   .. tab-item:: Python MPI example

      .. code-block:: python

         from mpi4py import MPI

         # Init MPI assuming that processes have been spawned
         comm = MPI.Comm.Get_parent()
         mpisize, mpirank, mpiroot = comm.Get_size(), comm.Get_rank(), 0
         if(mpirank == mpiroot) {
            # Read data from disk
         }
         # Copy inputs from the root to other processes
         data = MPI.COMM_WORLD.bcast(data, mpiroot=0)
         # Init inputs and simulate
         ...
         if(mpirank == mpiroot) {
            # Post-process data
         }

   .. tab-item:: C MPI example (for developers)

      .. code-block:: C

         #ifdef MPI
         #include <mpi.h>
         /* Init MPI */
         MPI_Init_thread(...)
         MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank);
         MPI_Comm_size(MPI_COMM_WORLD, mpi_size);

         /* Gather data from other processes to the root */
         if(mpi_rank == mpi_root) {
            MPI_Recv(data, datasize, datatype, id_sender, ...);
         else {
            MPI_Send(data, datasize, datatype, mpi_root, ...);
         }
         #endif

Threads (CPU)
*************

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

In some cases, sharing a node with multiple MPI processes can result in faster execution.
This could be due to how the cores access the memory or how threads are pinned to cores.
It is advisable to experiment with different configurations by adjusting the number of processes while keeping the total number of cores (``n x c``) constant.
This can help determine the optimal setup for your specific simulation.

*For Developers:* Most of the code runs within an OpenMP parallel region, with the primary exception being the input initialization.
It is crucial to ensure that the code remains thread-safe to prevent race conditions.
The example below demonstrates how to use OpenMP to spawn threads and update a shared variable safely, ensuring that only one thread modifies the variable at a time:

.. code-block:: C

   #include <omp.h>

   #pragma omp parallel
   {
      /* The brackets are optional if only one statement follows */
      dostuff(...);

      /* Ensure that only a single thread updates the variable at a time */
      #pragma omp atomic
      shared_variable += 1;
   }

Vectorization (CPU)
*******************

During the simulation, a group of markers belonging to the same thread is represented using a struct of arrays, where each array corresponds to a coordinate and the index within the array identifies which marker it belongs to.
The purpose of this struct is to locate markers in contiguous memory, enabling efficient vector operations on them.
These vector operations are executed within SIMD (Single Instruction Multiple Data) loops.

The number of markers within a vector is specified by the compiler flag NSIMD, which should be set to 8 for modern processors that support avx512 instructions.
Note that setting NSIMD=1 effectively disables vectorization.
It is also important to consider that vectorization introduces overhead; thus, one should not expect a linear increase in simulation efficiency with a higher NSIMD value, unlike the behavior observed with threads or processes.

The number of markers within a vector is specified by the compiler flag ``NSIMD``, which should be set to 8 for modern processors that support ``avx512`` instructions.
Note that setting ``NSIMD=1`` effectively disables vectorization.
It is also important to consider that vectorization introduces overhead; thus, one should not expect a linear increase in simulation efficiency with a higher ``NSIMD`` value, unlike the behavior observed with threads or processes.

Vectorization is implemented using OpenMP.
However, due to the complexity of the loops involved, only the Intel compiler `icc` has been found to reliably vectorize the code.
Typically, vectorization is enabled when the ``-march=<targetcpu>`` is set, provided that the target CPU supports vector instructions.
Additionally, you may need to set the flags ``-vecabi=cmdtarget`` or ``-ipo``.

*For Developers:* The basis of the vectorization lies in the structs of arrays, as demonstrated below:

.. code-block:: C

   typedef struct {
      real r[NSIMD];       /**< Particle R coordinate [m]          */
      real phi[NSIMD];     /**< Particle phi coordinate [phi]      */
      real z[NSIMD];       /**< Particle z coordinate [m]          */
      ...
   } particle_simd_fo;

Note that these arrays are dynamically allocated in the code since the same structs are also used in the GPU variant.
The SIMD loops are declared using the macro ``GPU_PARALLEL_LOOP_ALL_LEVELS``, which expands to ``#pragma omp simd`` when compiled for the CPU:

.. code-block:: C

   GPU_PARALLEL_LOOP_ALL_LEVELS
   for(int i = 0; i < NSIMD; i++) {
      if(p->running[i]) {
         p->z[i] += increment[i];
         dostuff(p->r[i], p->phi[i], p->z[i], &stuff);
      }
   }

To achieve effective vectorization, the code within the loop should not contain branches.
This means that:

- All ``for`` and ``while`` loops, as well as ``switch-case`` statements, must behave identically for all markers.
- There should be no multiple return statements in any functions that are called within the loop.

Additionally, you must instruct the compiler to generate vectorized versions of any functions called within a SIMD block.
This is done in the function declaration, as shown below:

.. code-block:: C

   GPU_DECLARE_TARGET_SIMD // For CPU expands to #pragma omp declare simd
   void dostuff(...);
   GPU_END_DECLARE_TARGET_SIMD // For CPU expands to nothing

Threads and Offloading (GPU)
****************************

In the GPU context, threads function differently from those on the CPU; they are not as independent.
Instead, GPU threads share similarities with vectorized operations, albeit with more flexibility.

GPU simulations always begin and end on the host CPU, which also manages all I/O operations.
Within ASCOT5, only the computationally intensive tasks that can be parallelized are executed on the GPU.
Since the CPU and GPU do not share memory, the simulation data must be offloaded from the host CPU to the target GPU.

To compile ASCOT5 for the GPU, set ``GPU=1``.
The GPU variant can be compiled using either OpenMP or OpenACC, so you must also set either ``OMP=1`` or ``ACC=1`` based on compiler support and effectiveness.

Although GPUs allow for the parallel simulation of many more markers than CPUs, computational performance may vary.
Whether the CPU or GPU is faster depends on the available hardware resources and the specifics of the simulation case.
Generally, only large simulations that fully utilize the GPU, i.e. the number of markers is hundreds of thousands at the very least, should be run on it when CPU nodes are also available.

*For Developers:* Since the CPU variant already supports vectorization and is thread-safe, the same code can be compiled for GPUs with minimal modifications.
These modifications are handled using custom macros that expand to different pragmas depending on whether the code is compiled for the CPU, GPU using OpenMP, or GPU using OpenACC.
To declare GPU parallel loops, use `GPU_PARALLEL_LOOP_ALL_LEVELS``, which is also used for SIMD loops in the CPU variant.

A key aspect specific to GPUs is offloading, where data is transferred between the host CPU and the GPU.
The following macros are used for offloading, expanding to the appropriate OpenMP or OpenACC pragmas:

.. code-block:: C

   /* Copy var from the host to the accelerator. Expands to:
    * - #pragma omp target enter data map(to:var) // OpenMP
    * - #pragma acc enter data copyin(var)        // OpenACC */
   GPU_MAP_TO_DEVICE(var)

   /* Copy var from the accelerator to the host. Expands to:
    * - #pragma omp target update from(to:var) // OpenMP
    * - #pragma acc update host(var)           // OpenACC */
   GPU_UPDATE_FROM_DEVICE(var)

In the code, inputs, diagnostics, and markers are offloaded before entering the simulation loop.
This implies that these variables must not be modified outside of GPU parallel loops, as modifications made within the GPU are not reflected on the CPU (and vice versa).
Only when the simulation ends are diagnostics and markers transferred back to the host CPU.

Finally, the compiler must be instructed to generate GPU-compatible versions of any function used within a GPU parallel loop.
This is done with the following macros:

.. code-block:: C

   /* Expands to:
    * - #pragma omp declare target // OpenMP
    * - #pragma acc routine seq    // OpenACC */
   GPU_DECLARE_TARGET_SIMD
   void dostuff(...);
   /* Expands to:
    * - #pragma omp end declare target // OpenMP
    * - nothing                        // OpenACC */
   GPU_END_DECLARE_TARGET

Accelerating simulations
************************

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
