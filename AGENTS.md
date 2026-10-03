# ASCOT5 Agent Notes

## Build and run

- The C solver is in `src/`; the Python package is `a5py/`. Root `make` moves compiled binaries/libraries to `build/`.
- Standard setup/build: activate `ascot-env`, then `make libascot -j`, `make ascot5_main -j`, and `pip install -e .`.
- Leonardo build from `~/ascot5`: load `nvhpc/25.11 hpcx-mpi/2.25.1 cuda/12.6 cmake/3.27.9 python/3.11.7 hdf5/1.14.3--hpcx-mpi--2.25.1--nvhpc--25.11`, then `srun -A fusio_hlst_6 -p boost_fua_prod --nodes=1 --ntasks=1 --cpus-per-task=16 --gpus-per-task=1 --time=00:10:00 make ACC=1 GPU=1 MPI=1 ascot5_main -j 16`. The root target puts the executable at `build/ascot5_main`.
- Leonardo batch runtime uses `$HOME/TorchFort/build/lib` in `LD_LIBRARY_PATH` and `$HOME/libatomic.so.1.2.0` in `LD_PRELOAD`; `launch_ascot.sh` loads those plus the matching HDF5 module.
- To run the training case after building: `cd /pitagora_scratch/userexternal/mpeybern/test_ascot5 && sbatch launch_ascot.sh`. For the validated warm start, use `BFLD_INITIAL_WEIGHTS=ft_div_1e-1_2075512.pt sbatch launch_ascot.sh`; without this variable the script initializes a fresh model. The script copies the config, creates `bfield_B_periodic_dB_init.pt`, and runs `ascot_fo_light` to `toto_dB_replay.h5`.
- Monitor with `squeue -j <jobid>` and `tail -f ascot-train-dB-<jobid>.out`; check completion with `sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed`.

## Training-model handoff

- Active work is on `feature/autograd_derivative_training`; integration commit is `2677ce1d` (on top of helper commits `c93d79d2` and `17688ebf`).
- Model source: `src/bfield_autodiff_model.py`; TorchFort config: `config_bfield_B_dB.yaml`; C integration and labels: `src/simulate/simulate_fo_fixed.c`, `src/bfield_nn.c`, `src/bfield_nn.h`; checkpoint loader/physical-unit inference: `src/wrap_bfield_autodiff.py`.
- Input is `x=((R-6)/2, (phi_mod-pi)/pi, Z/2)`. The model internally maps phi to `(sin(18 phi), cos(18 phi))`; ASCOT B_3DS confirms an 18-fold, 20-degree period. Keep C preprocessing and Python model normalization synchronized.
- B output scales are `(1.953079, 17.53601, 3.749450) T`. C and Python also share nine gradient-normalization constants and divergence scale `3.03743887`; if changing these, update both sides and the inference wrapper.
- TorchFort uses a trainable TorchScript SiLU model because its native MLP hardcodes ReLU. The model returns 13 targets: three normalized B components, nine autograd gradients (supervised by ASCOT `B_dB`), and weighted cylindrical divergence. Gradient forward uses `create_graph=True` so the derivative loss updates model weights.
- Each fixed-step simulation iteration trains with up to 2048 points: half current marker states, half deterministic uniform R-Z-phi replay points labelled by native `B_field_eval_B_dB`. The model is currently auxiliary training/validation only; particle pushing still uses native ASCOT fields. Do not describe it as having replaced the interpolation.
- The validated full run was Slurm job `2088692` (exit 0, 2m48s). On 10k independent uniform points, integrated-model B NRMSE was `(0.0576, 0.00717, 0.0467)` for `(B_R,B_phi,B_Z)`; the nine derivative NRMSEs were `(0.197,0.335,0.219,0.065,0.245,0.401,0.170,0.252,0.180)`. Its 20-degree periodicity error was below `7.7e-6 T`; normalized divergence RMS was `0.0358` vs native `8.5e-5`. `B_R` accuracy and divergence remain tradeoffs; boundary test B_R NRMSE was `0.212`.
- Detailed equations/results: `/pitagora_scratch/userexternal/mpeybern/test_ascot5/bfield_property_test_report_2075512.tex`; machine-readable validation: `integrated_dB_replay_validation_2088692.npz`; full run log: `ascot-train-dB-2088692.out`.
- Existing untracked items include `ascot-env/`, downloaded files, generated `.pt` artifacts, an OpenCode session JSON, and editor lock files. Do not stage or clean them wholesale; stage only intentional source/config changes.

## Python tests

- Unit/physics suites are `a5py/testascot/unittests.py` and `a5py/testascot/physicstests.py`; run them from `a5py/testascot/`. Unit tests need compiled binaries under `build/`.
