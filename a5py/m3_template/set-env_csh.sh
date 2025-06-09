
module purge
module load cineca



module load imasenv/3.42.0/foss/2023b/rc
#module unload HDF5/1.14.3-gompi-2023b
#module load HDF5/1.14.4.3-gompi-2023b
module load IMASPy/1.2.0-foss-2023b
module load mpi4py/3.1.5-gompi-2023b

module switch MUSCLE3/0.8.0-foss-2023b
module switch HDF5/1.14.4.3-gompi-2023b

setenv TEST_DIR $PWD
setenv PYTHONPATH /gss_efgw_work/work/g2diy/public/ascot/test_branch\:$PYTHONPATH
setenv PYTHONPATH /gss_efgw_work/work/g2diy/public/ascot/test_branch/venv/lib/python3.11/site-packages\:$PYTHONPATH
