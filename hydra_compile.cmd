#!/bin/bash -x

# @ shell=/bin/bash
# @ error = job.err.$(jobid)
# @ output = job.out.$(jobid)
# @ job_type = parallel
# @ requirements = (Feature=="mic")
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 1
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 00:11:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue

# load runtime libraries necessary to offload code to the Xeon Phi device
#module load intel/13.1 mkl/11.0 xeon-phi-rt/13.1

export OMP_NUM_THREADS=1
export KMP_AFFINITY=verbose,compact
export OFFLOAD_DEVICES=0
export OFFLOAD_REPORT=2
export MIC_ENV_PREFIX=MIC_
export MIC_OMP_NUM_THREADS=10
export MIC_KMP_AFFINITY=explicit,granularity=fine,proclist=[1-236:1]

# run the program

#cd /ptmp/${USER}/
#poe /u/${USER}/mictest.exe

#module unload intel
#module unload intel/16.0
module load   intel/16.0  hdf5-serial/1.8.16 xeon-phi-rt/16.0 gcc
module list

#make clean

#make TARGET=1 ascot5_gc FLAGS=" -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -fPIC -L/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib /u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib/libhdf5_hl.a /u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib/libhdf5.a -lrt -lz -ldl -lm -Wl,-rpath -Wl,/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib -I/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/include/"

#make clean
#make TARGET=1 hosts

#make TARGET=1 ascot5_gc FLAGS=" -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -fPIC -L/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib -lhdf5_hl -lhdf5 -lrt -lz -ldl -lm -Wl,-rpath -Wl,/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib -I/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/include/"

#nm /u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib/libhdf5.a > libhdf5.a.txt

make clean
#h5cc -c -o ascot4_interface.o ascot4_interface.c -fopenmp -qno-offload -std=c99 -Wall -DCOULOMBCOLL=0

#make TARGET=1 NSIMD=16 test_offload

#echo ----------------

#./test_offload

#echo ----------------



make TARGET=1 NSIMD=16 CC=h5cc ascot5_gc 

#make TARGET=1 ascot4_interface.o
#make TARGET=1 NSIMD=16 objs

#icc -o ascot5_gc ascot5_gc.o ascot4_interface.o hdf5_helpers.o hdf5_histogram.o hdf5_particlestate.o B_GS.o math.o consts.o wall_2d.o distributions.o B_2D.o B_ST.o B_TC.o plasma_1d.o interact.o step_gc_rk4.o step_fo_lf.o step_fo_vpa.o B_3D.o simulate_fo_lf.o simulate_gc_rk4.o wall_3d.o list.o octree.o particle.o endcond.o B_field.o E_field.o wall.o simulate.o orbit_write.o step_gc_cashkarp.o phys_orbit.o -fopenmp -std=c99 -Wall -DNSIMD=16 -DCOULOMBCOLL=0  -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -fPIC -L/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib /u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib/libhdf5_hl.a /u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib/libhdf5.a -lrt -lz -ldl -lm -Wl,-rpath -Wl,/u/system/SLES11/soft/hdf5/1.8.16/intel16.0/serial/lib
