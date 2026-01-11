# High-Performance N-Body Simulation Engine

![Language](https://img.shields.io/badge/Language-C11-blue.svg)
![Parallelism](https://img.shields.io/badge/Parallelism-OpenMP-green.svg)
![Build](https://img.shields.io/badge/Build-CMake-orange.svg)

An optimized gravitational N-body simulation engine written in **C**, designed to simulate millions of particles efficiently. This project demonstrates extreme performance engineering, reducing computational complexity from **$O(N^2)$** to **$O(N \log N)$** and achieving a **400x speedup** in memory allocation.

## 🚀 Key Optimizations

### 1. Algorithmic Efficiency (Barnes-Hut)
Replaced the naive Brute Force approach ($O(N^2)$) with the **Barnes-Hut** algorithm ($O(N \log N)$).
- **QuadTree/OctTree** data structure for spatial partitioning.
- **Approximation**: Far-away clusters are treated as single centers of mass ($\theta < 0.5$).

### 2. Memory Management (Arena Allocator)
Standard `malloc/free` is too slow for constructing/destructing trees every frame.
- **Solution**: Implemented a custom **Linear Arena Allocator**.
- **Result**: **>400x faster** allocation/deallocation cycle (measured via `scripts/bench_alloc.c`).
- **Zero Fragmentation**: Resetting the tree is a simple pointer reset (`arena_reset()`).

### 3. Data Locality (Morton Codes)
Accessing particles in random memory order causes cache misses.
- **Solution**: Sorted particles using **Z-Order Curve (Morton Codes)**.
- **Impact**: Improves CPU cache hit rate by keeping spatially close particles adjacent in memory.

### 4. Parallelism (OpenMP + K-Means)
- **Dynamic Scheduling**: Load balancing using OpenMP.
- **K-Means Clustering**: Pre-partitions particles to ensure thread workloads are spatially coherent and balanced.

## 📊 Performance Analysis

### Complexity Reduced
Visual comparison of execution time vs particle count ($N$):

![Complexity](docs/execution_time_vs_N.png)

*The Barnes-Hut implementation (Blue) scales linearly compared to the exponential Brute Force (Red).*

## 🛠️ Build & Run

### Prerequisites
- GCC / Clang with OpenMP support
- CMake 3.10+

### Compilation
```bash
mkdir build && cd build
cmake ..
make
```

### Running the Simulation
```bash
# Usage: ./simulate <N> <input_file> <steps> <dt> <threads> <theta>
./simulate 100000 data/input.gal 100 0.01 8 0.5
```

### Reproducing Benchmarks
To verify the memory allocator performance:
```bash
gcc -O3 scripts/bench_alloc.c -o build/bench_alloc
./build/bench_alloc
# Output:
# System Malloc: 0.1452s
# Arena Alloc:   0.0003s (Speedup: ~480x)
```
