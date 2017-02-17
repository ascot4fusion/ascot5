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
export MIC_OMP_NUM_THREADS=1
export MIC_KMP_AFFINITY=explicit,granularity=fine,proclist=[1-236:1]

# run the program

#cd /ptmp/${USER}/
#poe /u/${USER}/mictest.exe

#module unload intel
module load   hdf5-serial/1.8.16 xeon-phi-rt/16.0 gcc
module list

make clean
make TARGET=1
