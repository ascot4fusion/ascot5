===============
Parallelization
===============

.. figure:: ../figures/ascot5parallel.svg
   :class: with-border
   :align: center
   :width: 300px

   Illustration of all the levels at which ASCOT5 is parallelized.

It is important to have a grasp on how the parallelization in ASCOT5 works, in order to fully utilize the available HPC resources.

ASCOT5 is a hybrid code that is parallelized using both MPI and OpenMP.
MPI (Message Passing Interface) distributes the simulation among different computing nodes while OpenMP (Open Multi-Processing) ensures all cores within each node are fully utilized.
The difference between MPI and OpenMP is that the memory is not shared between different *MPI processes* whereas all OpenMP "processes" (*threads*) access the same memory.

The figure above illustrates how the code is parallelized.
Let's go through these points in detail:

1. At the beginning when the simulation is launched, `n` identical processes are started with each running ASCOT5.

2. Each MPI process reads the data from the disk and prepares it for the simulation (interpolation splines are constructed, wall octree is formed, etc.).

3. Each MPI process reads the marker data into a single (huge) queue.
   Up to this point all MPI processes have performed identical tasks.

4. Now each MPI process picks a chunk of the marker queue.
   For example, assuming we have 2999 markers and three MPI processes, the root process picks first 1000, the second process picks markers 1001-2000, and the third and final process picks the remaining 999.

5. Each MPI process initializes the markers it has picked. Then the root process gathers all markers from all MPI processes (temporarily) to write the inistate on disk.

6. Now *offloading* happens.
   More about that later.

7. Everything is ready for the simulation and all that happens now is identical between the MPI processes with the exception that each process has its own marker population to simulate.

8. OpenMP spawns threads that is equal to the amount of CPUs on the node.
   Since data is shared between the threads, every thread can interpolate the same magnetic field data and update the same distribution.
   This is why it is preferred to have a single MPI process per node because then one can use all available memory for collecting high-resolution distributions.
   If two or more MPI processes were run in the same node, then each would reserve space for that processes' input data and distribution since no data is shared between MPI processes while markers are being simulated.

9. Each thread picks markers from the marker queue that the MPI process (where those threads belong to) is supposed to simulate.
   Every time an individual thread finishes simulating a marker, it picks a new marker from the queue until none is left.

10. But wait!
    That's not all.
    One thread does not simulate just a one marker at the time.
    Instead, it simulates a vector of markers or *SIMD array* of markers.

11. SIMD (Single Instruction, Multiple Data) array differs from threads in that while each thread works independently, in the SIMD array all operations are performed at the same time for all the markers in the array.
    For example, if there are eight markers in the array, then collisions are evaluated concurrently, orbit integration is performed concurrently etc for those eight markers.
    When a marker in an array has finished simulation, the thread picks a new marker from the queue to replace it.

12. So to summarize: markers are first distributed evenly among the MPI processes which are independent processes that share no memory.
    Then within each MPI process threads are spawned and each thread works independently, picking new markers from the queue every time that thread has finished simulating one.
    Threads share the input and output data (memory) but they work in their own CPUs.
    Modern CPUs can perform vector operations, meaning that each thread can have a SIMD array of markers which are pushed concurrently.
    So in total we have three layers on which the code is parallelized.

13. Once each thread has finished and none of the MPI processes have unfinished markers, the root process collects the simulated markers and all diagnostics from all the MPI processes, and writes those to a disk.

Offloading
==========
