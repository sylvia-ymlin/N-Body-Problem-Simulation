# High‑Performance N‑Body Simulation Engine

An optimized gravitational N‑body simulation engine written in C, leveraging the **Barnes–Hut** algorithm to reduce complexity from $O(N^2)$ to $O(N \log N)$. This project focuses on High-Performance Computing (HPC) techniques including custom memory allocators, spatial sorting, and multicore parallelism.

## Project Structure

```
├── src/                # Source Code
│   ├── galsim.c        # Main simulation loop & Arg parsing
│   ├── barnes_hut.c    # QuadTree construction & Force calculation (Core)
│   ├── kmeans.c        # K-means clustering for load balancing
│   ├── morton.c        # Z-Order (Morton Code) sorting for cache locality
│   ├── ds.h            # Data structures & Arena definition
│   └── legacy/         # Archived original implementations (v1)
├── docs/               # Documentation
│   ├── interview_prep.md    # Technical deep-dive & optimization story
│   └── PROJECT_ANALYSIS.md  # Theoretical analysis
├── scripts/            # Helper Scripts
│   ├── generate_data.py     # Generate binary input files
│   ├── bench_alloc.c        # Micro-benchmark for Memory Arena
│   └── bench_traversal.c    # Micro-benchmark for Tree Traversal
├── data/               # Input/Output data
├── CMakeLists.txt      # Build configuration
└── README.md           # This file
```

## Key Optimizations

1.  **Algorithmic**: Barnes-Hut Tree Code ($O(N \log N)$).
2.  **Memory Management**: Custom **Linear Arena Allocator** replaces `malloc/free` for tree nodes, achieving >400x speedup in allocation.
3.  **Data Locality**: 
    *   **Z-Order Sorting (Morton Code)**: Reorders particles in memory to match geometric proximity (Cache Friendly).
    *   **K-Means Clustering**: Partitions particles for balanced OpenMP scheduling.
4.  **Parallelism**: OpenMP Dynamic Scheduling for force calculation and K-means.
5.  **Compute Kernel**: **Stackless Tree Traversal** implementation reduces function call overhead.
6.  **Safety**: Removed all Variable Length Arrays (VLAs) to support simulation of millions of particles without Stack Overflow.

## Build

Use CMake to build the project (requires OpenMP):

```bash
mkdir -p build
cd build
cmake ..
make
```

## Usage

### Run Simulation
```bash
# Syntax
./simulate <N> <filename> <nsteps> <delta_t> <n_threads> <theta_max> <k_clusters>

# Example: 100k particles, 100 steps, using 4 threads and 10 clusters
./simulate 100000 ../data/input.gal 100 0.001 4 0.5 10
```

### Run Benchmarks
You can compile and run specific micro-benchmarks located in `scripts/`:

```bash
# Benchmark Memory Arena vs System Malloc
gcc -O3 scripts/bench_alloc.c -o build/bench_alloc
./build/bench_alloc

# Benchmark Recursive vs Stackless Traversal (requires linking core objs)
# (See scripts/bench_traversal.c for details)
```

## Input File Format
Binary file containing $N$ particles sequentially. Each particle consists of 6 `double` values:
`[pos_x, pos_y, mass, vx, vy, brightness]`
