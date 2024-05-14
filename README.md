# N-Body Simulation in C

This project simulates a 2D gravitational N-body system and compares several increasingly faster force-computation kernels. The code starts from a simple all-pairs baseline and moves to a Barnes-Hut based implementation with memory-layout and OpenMP improvements.

The project is intentionally small in scope. Its main purpose is to show a clear optimization path for one problem, not to be a general-purpose physics engine.

## What Problem It Solves

Direct N-body simulation is expensive because every particle interacts with every other particle. A naive implementation costs `O(N^2)` work per force update, which becomes impractical as `N` grows.

This repository shows a compact progression from:

- exact pairwise force calculation
- to Barnes-Hut approximation
- to better memory locality
- to shared-memory parallelism with OpenMP

## General Approach

The executable supports five versions:

- `v1`: naive all-pairs reference implementation
- `v2`: Barnes-Hut tree to reduce force computation cost
- `v3`: Barnes-Hut with arena allocation for tree nodes
- `v4`: Barnes-Hut with Morton ordering for better locality
- `v5`: `v4` plus OpenMP parallel force evaluation

The main workflow is the same for every version:

1. read particle data
2. compute forces with the selected kernel
3. advance the system with Velocity Verlet integration
4. write the final particle state to `data/outputs/`

`v4` also contains an optional K-means reordering path. It is kept as an experiment for comparison, but it is not the main recommended configuration.

## Project Structure

```text
.
├── main.c                  # command-line entry point and simulation loop
├── simulation.h            # shared kernel interface
├── v1_naive.c              # O(N^2) reference kernel
├── v2_barnes_hut.c         # Barnes-Hut implementation
├── v3_arena.c              # Barnes-Hut with arena allocation
├── v4_morton.c             # Barnes-Hut with Morton ordering
├── v5_parallel.c           # OpenMP version of v4
├── core/
│   ├── io.c/.h             # binary input/output
│   ├── domain.c/.h         # simulation bounding box helpers
│   ├── morton.c/.h         # Z-order sorting
│   ├── kmeans.c/.h         # optional clustering experiment
│   ├── ds.h                # tree and arena data structures
│   ├── types.h             # shared structs
│   └── time_utils.h        # timing helper
├── scripts/
│   ├── generate_data.py    # create `.gal` input files
│   ├── benchmark_suite.py  # strong/weak scaling runs
│   ├── benchmark_full_suite.py
│   ├── benchmark_ablation.py
│   ├── accuracy.py         # compare result files
│   └── plot_scaling.py     # create figures from metrics
├── data/
│   ├── inputs/             # generated inputs
│   ├── outputs/            # simulation outputs
│   └── metrics/            # saved benchmark results
└── figures/                # generated plots used in the report
```

## Build

Requirements:

- C compiler with C99 support
- CMake
- OpenMP support if you want to run `v5`
- Python 3 for data generation and benchmark scripts

Build the project:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

Generate an input file:

```bash
python3 scripts/generate_data.py 100000 data/inputs/input_100k.gal disk
```

Run the simulator:

```bash
./build/nbody_simulate 5 100000 data/inputs/input_100k.gal 1 1e-5 8 0.5 0
```

Arguments:

```text
./build/nbody_simulate <version> <N> <input.gal> <nsteps> <dt> <n_threads> <theta> <k>
```

- `version`: kernel version `1` to `5`
- `N`: number of particles to read
- `input.gal`: binary input file
- `nsteps`: number of integration steps
- `dt`: timestep size
- `n_threads`: OpenMP thread count
- `theta`: Barnes-Hut acceptance threshold
- `k`: number of K-means clusters for `v4`; use `0` for normal Morton ordering

## Example Output

A typical run prints the selected version, progress, and total runtime. `v5` also reports tree-build and force-computation timings.

```text
Running v5 | N=100000 | Steps=1 | Threads=8
[HP] Schedule: dynamic,128 | Tree Build: 0.07 s | Force Compute: 0.05 s
Step 0/1
Simulation Complete: 0.17s
```

The result is written to `data/outputs/result_v5.gal`.

## Benchmark Scripts

The benchmark scripts are optional. They are useful for reproducing the saved metrics and figures, but they are not required to understand or run the simulator itself.

```bash
python3 scripts/benchmark_full_suite.py
python3 scripts/benchmark_ablation.py
python3 scripts/plot_scaling.py
```

The repository already includes example metric files in `data/metrics/` and plots in `figures/`.

## Results Summary

The repository includes saved benchmark data in `data/metrics/` and plots in `figures/`. The main results are:

### Serial Versions

The biggest improvement comes from switching from exact all-pairs force calculation to Barnes-Hut. Arena allocation helps a little, while Morton ordering gives the largest implementation-level gain after the algorithm change.

| Version | Main idea | Example runtime |
| :------ | :-------- | --------------: |
| `v1` | naive all-pairs | `86.11s` |
| `v2` | Barnes-Hut | `11.06s` |
| `v3` | Barnes-Hut + arena | `10.66s` |
| `v4` | Barnes-Hut + Morton | `5.13s` |

These numbers are from the included serial comparison at `N=10k`, `nsteps=200`, `dt=1e-5`.

### Parallel Version

For `v5`, the force loop scales well at first, but total speedup eventually levels off because tree construction and Morton ordering are still serial.

At `N=100000`, the saved strong-scaling run reaches about `5.1x` end-to-end speedup by `16` to `20` threads. After that, the serial part of the program becomes the main bottleneck.

### Correctness Check

The project also includes a simple accuracy check against the naive reference. When Barnes-Hut is run with `theta=0`, it matches the reference solution up to floating-point rounding error. With the default `theta=0.5`, the saved error remains small enough for this project’s engineering use case.

## Limitations

- The simulator is 2D only.
- The code focuses on one optimization path rather than a broad set of features.
- Tree construction is still serial, so parallel speedup is limited.
- Input and output use a project-specific binary format.
- The K-means path is experimental and mainly useful for comparison, not as the default method.
- This is a study project, not a production-grade simulation package.
