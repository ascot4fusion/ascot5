# ASCOT5 Agent Notes

- The C simulation is in `src/`; the Python package is `a5py/`. The root Makefile builds C targets and moves executables/libraries into `build/`.
- Standard setup/build: `conda env create -f environment.yaml`, activate `ascot-env`, then `make libascot -j`, `make ascot5_main -j`, and `pip install -e .`.
- On Leonardo with the training branch, load `nvhpc/25.11 hpcx-mpi/2.25.1 cuda/12.6 cmake/3.27.9 python/3.11.7`; also load `hdf5/1.14.3--hpcx-mpi--2.25.1--nvhpc--25.11` at runtime. The executable needs `$HOME/TorchFort/build/lib` in `LD_LIBRARY_PATH` and `$HOME/libatomic.so.1.2.0` in `LD_PRELOAD`.
- Leonardo training build from `src/`: `srun make ACC=1 GPU=1 MPI=1 ascot5_main -j 16`.
- The `feature/training_model` path uses `config_bfield_B.yaml` and trains B only. Derivative inference is computed from a saved B model with autograd; it is not derivative-target training. The wrapper returns derivatives with respect to normalized coordinates/field outputs, so convert scales before comparing with ASCOT physical derivatives.
- User simulation tests run in `/pitagora_scratch/userexternal/mpeybern/test_ascot5`; run `sbatch launch_ascot.sh` there. The script expects the source executable at `~/ascot5/src/ascot5_main` and uses `ascot_fo_light` as input.
- Python unit/physics suites are `a5py/testascot/unittests.py` and `a5py/testascot/physicstests.py`; the unit suite expects compiled test binaries under `build/`.
