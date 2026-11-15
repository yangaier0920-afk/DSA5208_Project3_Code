# DSA5208 Project 3: Jacobi solvers for the 2D/3D Poisson equation

This package contains the source code of the project and the results shown in the
report. All programs solve the manufactured Poisson problem −Δu = f on the unit
square/cube with Dirichlet boundary values, using Jacobi sweeps on the 5-point (2D)
or 7-point (3D) stencil: a serial solver, a Cartesian MPI and MPI+OpenMP solver, and
an OpenMP CPU / OpenMP target (GPU) solver.

```text
README.md                     this file
DSA5208_Project3_Report.pdf   the project report
requirements.txt              Python packages for the analysis and plotting scripts
code/                         source code, runners, validation and analysis scripts, PBS job scripts
  Makefile                    builds all CPU programs
  serial/                     serial solver: convergence/accuracy study and problem-size scaling
  mpi/                        Cartesian MPI and MPI+OpenMP solver, scaling runs, diagnostics
  openmp/                     OpenMP CPU and OpenMP target (GPU) solver
output/                       results shown in the report
  figures/                    the ten result figures (fig01 ... fig10)
  tables/                     summary tables behind the figures and tables of the report
```

## Code

| Path (under `code/`) | Purpose |
| --- | --- |
| `serial/jacobi_numerics.cpp` | Serial solver with fixed-work and convergence (solve) modes, used for the convergence and accuracy study |
| `serial/run_a_numerics.py`, `analyze_a_numerics.py`, `vanda_a_numerics.pbs` | Runner, analysis and PBS script for that study |
| `serial/scaling/` | Serial problem-size scaling: solver (`jacobi_a_scaling.cpp`), runner, analysis and PBS script |
| `serial/historical_2d/` | Historical 2D programs (`jacobi.cpp`, `jacobi_serial_local_update.cpp`, `jacobi_serial_row_update.cpp`, `jacobi_mpi.cpp`) with the driver and PBS script for the serial versus single-rank MPI comparison |
| `mpi/jacobi_cart.cpp` | Cartesian MPI solver; compiled with OpenMP it is the MPI+OpenMP (hybrid) solver |
| `mpi/run_b.py`, `exec_affinity.py`, `verify_b.py`, `analyze_b.py`, `plot_b.py`, `make_jobs.py` | Strong/weak-scaling runner, CPU-binding wrapper, validation, analysis and job generator |
| `mpi/jobs/` | PBS scripts for P = 1, 8, 12, 36, 72, 144, 288, 432, 720 and a multi-node validation job |
| `mpi/supplementary/` | Instrumented solver variants (`jacobi_cart_tuned.cpp`, derived from `jacobi_cart.cpp` by `prepare_solver.py`), two-rank ping-pong (`mpi_link_probe.cpp`), runner, variant validation and analysis |
| `mpi/supplementary/jobs/` | PBS scripts for the exclusive single-core baseline, the 144-core runs, the 72-core hybrid ablation, the ping-pong and the 2D MPI check |
| `openmp/jacobi_c_openmp.cpp`, `jacobi_c_target.c`, `jacobi_c_target.h` | OpenMP CPU solver and OpenMP target (GPU) kernel |
| `openmp/build.sh`, `run_c.py`, `make_job.py`, `analyze_c.py`, `cluster_env.sh.example` | Build script, runner, PBS job generator, analysis and NVHPC environment |
| `openmp/followup/` | Same-node single-thread comparison, GPU data-region phase timing, target-call probe and CPU thread scaling |
| `vanda_modules.example.sh`, `mpi/runtime_modules.sh` | Vanda module setup (GCC 12.2.0, Open MPI 4.1.5, Python 3.10) |

## Requirements

- Linux with a C++17 compiler (GCC 12.2.0 for the reported CPU runs) and Python 3.
  The runners and validation scripts use only the Python standard library; the
  analysis and plotting scripts need matplotlib (`pip install -r requirements.txt`).
- MPI with `MPI_THREAD_FUNNELED` support (Open MPI 4.1.5 for the reported runs).
- GPU build: NVIDIA HPC SDK 24.7 (CUDA 12.3) and an NVIDIA GPU (an A40 was used).
- On Vanda, load the modules first (from `code/`):
  `source mpi/runtime_modules.sh` for the CPU programs, or
  `cp openmp/cluster_env.sh.example openmp/cluster_env.sh && source openmp/cluster_env.sh` for the GPU program.
  Compile and run on a compute allocation, not on a login node.

## Build

All commands below are run from the `code/` directory.

```sh
cd code
make                                              # all CPU programs
make numerics scaling historical mpi openmp-cpu   # or individual groups
bash openmp/build.sh nvidia                       # GPU program (NVHPC, on a GPU node)
```

`make mpi` builds `build/jacobi_cart_mpi` (pure MPI), `build/jacobi_cart_hybrid`
(MPI+OpenMP), `build/jacobi_cart_tuned` and `build/jacobi_cart_axis` (instrumented
variants with and without OpenMP) and `build/mpi_link_probe`. `make openmp-cpu`
builds `openmp/build/jacobi_c_cpu`. The GPU build compiles the C++ host code and
the C target code separately and links them into `openmp/build/jacobi_c_gpu`.

## Correctness checks

Each check takes seconds and writes to a new output directory.

```sh
# Serial solver: 7x7x7 manufactured problem, error against the exact discrete solution ~1e-15
build/jacobi_numerics --problem poisson3d_mms --nx 7 --mode fixed --max-iters 1000
python3 serial/run_a_numerics.py --binary build/jacobi_numerics --suite smoke --out-dir check_numerics
python3 serial/scaling/run_a_scaling.py --binary build/jacobi_a_scaling --out-dir check_scaling --smoke

# MPI / hybrid solver: six oracle cases (uneven partitions, 2D and 3D, hybrid threads)
python3 mpi/verify_b.py --binary build/jacobi_cart_hybrid --launcher "mpirun --bind-to none" --out check_mpi
python3 mpi/supplementary/validate_variants.py --binary build/jacobi_cart_tuned --launcher "mpirun --bind-to none" --out check_variants

# OpenMP CPU solver: discrete-solution, transient and boundary checks
python3 openmp/run_c.py --binary openmp/build/jacobi_c_cpu --suite check --threads 2 --cpu-only --out check_openmp
```

The MPI checks start up to eight ranks. On a machine with fewer cores add
`--oversubscribe` to the launcher and set `OMP_WAIT_POLICY=PASSIVE`, for example
`--launcher "mpirun --oversubscribe --bind-to none -x OMP_WAIT_POLICY"`.
`openmp/run_c.py --suite check` without `--cpu-only`, run with the GPU binary on a
GPU node, adds the CPU/GPU full-grid comparisons and checks that the target region
runs on the device.

## Running the solvers directly

Every solver prints one JSON line with the timing, the residual and the errors
against the exact solutions. Timed loops contain only the sweeps (and, in parallel
runs, the halo exchange); initialization and final diagnostics are outside the timer.

```sh
# Serial (fixed work of K sweeps, or solve to a relative residual tolerance)
build/jacobi_a_scaling --problem poisson3d_mms --nx 128 --mode fixed --max-iters 1000
build/jacobi_numerics  --problem poisson3d_mms --nx 31 --mode solve --max-iters 100000 --tol 1e-10 --check-every 10

# MPI: process grid px*py*pz must equal the number of ranks
mpirun -np 8 build/jacobi_cart_mpi --dim 3 --nx 256 --ny 256 --nz 256 --px 2 --py 2 --pz 2 --iterations 1000

# MPI+OpenMP: 12 ranks x 6 threads
OMP_NUM_THREADS=6 mpirun -np 12 --bind-to none -x OMP_NUM_THREADS \
    build/jacobi_cart_hybrid --dim 3 --nx 256 --ny 256 --nz 256 --px 3 --py 2 --pz 2 --iterations 1000

# Instrumented variant with per-rank phase and per-axis timers
# (--variant baseline|parallel_pack|first_touch|persistent|optimized)
OMP_NUM_THREADS=6 mpirun -np 12 --bind-to none -x OMP_NUM_THREADS \
    build/jacobi_cart_tuned --variant parallel_pack --dim 3 --nx 256 --ny 256 --nz 256 \
    --px 3 --py 2 --pz 2 --iterations 1000 --diagnostic 1

# OpenMP CPU and GPU
openmp/build/jacobi_c_cpu --backend cpu --dim 3 --nx 256 --iterations 1000 --threads 12
OMP_TARGET_OFFLOAD=MANDATORY openmp/build/jacobi_c_gpu --backend gpu --dim 3 --nx 256 --iterations 1000 --threads 12

# Historical 2D program: interior size, sweep count, tolerance (0 = exactly that many sweeps)
build/jacobi_serial 1024 1000 0
```

For the 2D MPI problem use `--dim 2 --pz 1`. The `rank_timings` field of a
diagnostic run lists, for every rank, the total, update, pack/unpack,
communication/wait and swap times.

## Cluster experiments (Vanda, PBS)

Replace `YOUR_PROJECT` in the job scripts with a valid project, check the queue and
resource lines for your account, and submit from the directory shown (paths are
relative to `code/`). Each job compiles the code on its allocation, records the
allocated hosts and writes its results under `results/` in the submission directory
(the OpenMP follow-up jobs write to `openmp/followup/runs/`).

| Experiment | Submit |
| --- | --- |
| Serial problem-size scaling | `cd serial/scaling && qsub vanda_a_scaling.pbs` |
| Convergence and accuracy study | `cd serial && qsub vanda_a_numerics.pbs` |
| Serial versus single-rank MPI (historical 2D) | `cd serial/historical_2d && qsub vanda_single_rank.pbs` |
| MPI strong and weak scaling at P cores | `qsub mpi/jobs/b_<P>.pbs` (regenerate with `python3 mpi/make_jobs.py --out mpi/jobs`) |
| Multi-node MPI validation | `qsub mpi/jobs/b_validate.pbs` |
| Exclusive single-core baseline (MPI P=1 and serial on one bound core) | `qsub mpi/supplementary/jobs/baseline.pbs` |
| 144-core runs (nine and four hosts, hybrid layouts) | `qsub mpi/supplementary/jobs/clean144.pbs`, `qsub mpi/supplementary/jobs/paired144.pbs` |
| 72-core MPI+OpenMP ablation | `qsub mpi/supplementary/jobs/optimize72.pbs` |
| Two-rank ping-pong | `qsub mpi/supplementary/jobs/link_reference.pbs` (edit the two host names) |
| 2D MPI scaling (512×512) | `qsub mpi/supplementary/jobs/twod.pbs` |
| CPU/GPU matrix (CPU1, CPU12, GPU at 64³/128³/256³) | see below |
| Same-node single-thread comparison, GPU phases, thread scaling | `cd openmp/followup && qsub -P YOUR_PROJECT -v C_SUITE=cpu followup.pbs` (also `C_SUITE=gpu`, `C_SUITE=scaling`) |

The CPU/GPU job is generated in the `openmp` directory:

```sh
cd openmp
cp cluster_env.sh.example cluster_env.sh
python3 make_job.py --queue auto --project YOUR_PROJECT --select '1:ncpus=36:ngpus=1' \
    --threads 12 --suite formal --walltime 01:00:00 --confirmed-resources --output c_formal.pbs
qsub c_formal.pbs
```

Measurement protocol used by all runners: a fixed work of K = 1000 sweeps, one
warmup run and three measured fresh-start runs per configuration (five for the
single-core baselines), all valid samples kept, and medians with the observed
minimum and maximum reported. MPI times are the maximum over ranks.

The analysis scripts read the output directories written by the runners:

```sh
python3 serial/scaling/analyze_a_scaling.py --output <figure_dir> <scaling_output_dir>
python3 serial/analyze_a_numerics.py <numerics_output_dir>
python3 mpi/analyze_b.py --input <run_b_output_root> --out <analysis_dir>
python3 mpi/supplementary/analyze_feedback.py --input <run_feedback_output_root> --out <analysis_dir>
python3 openmp/analyze_c.py <openmp_output_dir>
python3 openmp/followup/plot_scaling.py <thread-scaling CSV from the scaling suite>
```

## Output

`output/` holds the results shown in the report: the ten result figures and the
summary tables they were drawn from. Re-running an experiment (the PBS jobs above,
followed by the analysis scripts) regenerates the per-run records, from which the
analysis scripts recompute these summaries and draw the corresponding plots (their
panel layout can differ from the report figures). The per-run records and job logs
themselves are not included in this package. The
tables listed below hold the numerical data behind every figure and the main result
tables of the report, so these values can be checked without re-running the
experiments. Times are in seconds (medians with the observed
minimum and maximum over the retained samples) unless the file states otherwise.
In the MPI tables, `run` is `first` for the first measurement campaign and
`supplementary` for the later runs on separate allocations.

| File | Content | Report |
| --- | --- | --- |
| `figures/fig01_serial_scaling.png` | Serial fixed-work time and cost per point update | Figure 1 |
| `figures/fig02_numerical_accuracy.png` | Residual histories and grid accuracy | Figure 2 |
| `figures/fig03_mpi_time.png` to `fig05_mpi_efficiency.png` | Strong-scaling time, speedup and efficiency | Figures 3–5 |
| `figures/fig06_weak_scaling.png` | Weak scaling (64³ points per rank) | Figure 6 |
| `figures/fig07_communication.png` | Phase shares and per-axis waits | Figure 7 |
| `figures/fig08_hybrid_layouts.png` | 72-core hybrid ablation and 144-core hybrid layouts | Figure 8 |
| `figures/fig09_gpu_runtime.png`, `fig10_gpu_speedup.png` | CPU/GPU runtimes and GPU speedups | Figures 9–10 |
| `tables/serial/serial_scaling.csv` | Serial scaling, 2D n=256–4096 and 3D n=16–256 | Section 3.1 |
| `tables/serial/numerical_runs.csv`, `grid_accuracy.csv`, `residual_history_n{15,31,63}.csv` | Convergence and accuracy study | Section 3.2 |
| `tables/mpi/single_core_baseline.csv` | MPI P=1 and serial program on one bound core (five samples) | Sections 3.3, 4.2 |
| `tables/mpi/strong_scaling.csv` | Strong scaling with speedup and efficiency | Section 4.3 |
| `tables/mpi/weak_scaling.csv` | Weak scaling with efficiencies against P=1, 8 and 72 | Section 4.4 |
| `tables/mpi/p144_allocations.csv` | The same P=144 workloads on different allocations, including the hybrid layouts | Sections 4.6, 5.1 |
| `tables/mpi/phase_breakdown_256.csv`, `p144_axis_waits.csv` | Phase breakdown of the slowest rank and per-axis Sendrecv waits | Sections 4.5–4.6, Appendix A |
| `tables/mpi/ping_pong.csv`, `twod_512_scaling.csv` | Two-rank round trips and the 2D MPI check | Section 4.5 |
| `tables/mpi/hybrid_ablation_72.csv`, `hybrid_layouts_144.csv` | 72-core ablation and 144-core hybrid layouts | Section 5.1 |
| `tables/openmp/cpu_gpu_formal.csv`, `cpu_gpu_speedups.csv` | CPU1/CPU12/GPU times and GPU speedups | Sections 5.3–5.4 |
| `tables/openmp/effective_throughput.csv` | Modeled effective throughput and data volumes | Sections 5.3–5.4 |
| `tables/openmp/same_node_single_thread.csv`, `thread_scaling.csv` | Same-node single-thread comparison and 1–36 thread scaling | Sections 3.3, 5.3 |
| `tables/openmp/gpu_phases.csv`, `target_call_probe.csv` | GPU data-region phases (seconds) and target-call probe | Section 5.4 |
