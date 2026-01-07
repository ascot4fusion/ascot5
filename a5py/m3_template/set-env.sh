module purge
module load cineca

module load imasenv/3.42.0/foss/2023b/rc

module load IMASPy/1.2.0-foss-2023b
module load mpi4py/3.1.5-gompi-2023b

module switch MUSCLE3/0.8.0-foss-2023b

export TEST_DIR=$( realpath $( dirname -- "${BASH_SOURCE[0]}" ) )
