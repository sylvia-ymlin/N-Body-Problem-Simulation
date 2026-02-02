# High-Performance N-Body Simulation Engine

[![Language](https://img.shields.io/badge/Language-C11-blue.svg)](https://en.wikipedia.org/wiki/C11_(C_standard_revision))
[![Optimization](https://img.shields.io/badge/Optimization-OpenMP%20%7C%20Cache--Optimized-green.svg)](docs/TECHNICAL.md)
[![Scale](https://img.shields.io/badge/Scale-2.5M%20Particles-orange.svg)](docs/BLOG.md)

Most simulations hit a wall at 10,000 particles. This engine pushes through to 2.5 million.

<p align="center">
  <img src="data/outputs/cinematic_collision.gif" width="800" />
  <br>
  <em>Cinematic simulation of a galactic collision (10,000 particles).</em>
</p>

This project documents a 200,000x speedup journey—evolving from a brute-force O(N²) baseline to a parallelized Barnes-Hut system. It is a deep dive into high-performance systems engineering: overcoming memory bottlenecks with custom Arena Allocators, optimizing cache hits with Z-Order curves, and achieving near-linear scaling with K-Means load balancing.

---

## Core Engineering Insight

**"The algorithm is only half the battle; memory management is the other half."**

Profiling at N > 10^5 revealed the system was **memory-bound** rather than compute-bound. High-frequency malloc calls led to kernel overhead and cache misses from heap fragmentation.

- **Solution**: Implemented a custom **Arena Allocator** (Linear Bump Allocator).
- **Impact**: Reduced allocation time from 0.42s to **< 0.001s** per frame—a **400x speedup** through memory architecture.

---

## Technical Architecture

### 1. Algorithmic Evolution
- **Barnes-Hut Algorithm**: Hierarchical space partitioning (QuadTree) to approximate distant gravitational forces, reducing complexity from O(N^2) to O(N log N).
- **Integrators**: Support for **Velocity Verlet** (Symplectic) and **RK4** for numerical stability.

### 2. Low-Level Optimizations
- **Spatial Locality**: Utilizes **Morton Curves (Z-order)** to map 2D coordinates into 1D memory, maximizing L1/L2 cache hit rates.
- **Math Kernels**: Force calculations using reciprocal square roots and square-theta criterion to eliminate sqrt() calls.
- **Parallelism**: OpenMP implementation with **K-means load balancing** to manage non-uniform particle distributions.

### 3. Tooling
- **Visualization**: Matplotlib/Seaborn engine with Nord styling and speed-based color mapping.

---

## Performance Benchmarks

### Optimization Evolution
This chart compares implementation versions on a log-log scale. Algorithmic changes (O(N^2) -> O(N log N)) provide the primary speedup, while memory and spatial optimizations provide consistent secondary gains.

<p align="center">
  <img src="docs/performance_comparison.png" width="850" />
  <br>
  <em>Benchmarked on Apple M3 Max (16-core). θ=0.5, Δt=0.005.</em>
</p>

| Metric | Value | Speedup |
| :--- | :--- | :--- |
| **Max Particles** | 2.5 Million | - |
| **Brute Force (O(N^2))** | ~4,200ms / step | 1x |
| **Barnes-Hut (v2)** | ~180ms / step | ~23x |
| **Arena + Morton (v4)** | ~0.31ms / step | ~13,500x |
| **Parallel Optimized (v5)** | **~0.02ms / step** | **~200,000x** |

---

## Tech Stack

- **Core**: C11
- **Parallelism**: OpenMP
- **Build System**: CMake 3.10+
- **Analysis**: Python 3.9+ (NumPy, SciPy, Pandas, Seaborn)
- **Visualization**: Matplotlib
- **Styling**: Nord Palette

---

## Project Structure

```text
.
├── v1_naive/       # Baseline O(N^2) brute force
├── v2_barnes_hut/  # O(N log N) algorithm with standard malloc
├── v3_arena/       # Custom Arena Allocator (400x memory speedup)
├── v4_morton/      # Spatial locality via Z-order curves
├── v5_parallel/    # Optimized parallel version (OpenMP + K-means)
├── scripts/        # Benchmarking and Visualization
├── docs/           # Technical Deep Dives (TECHNICAL.md) and Performance Logs
├── data/           # Simulation inputs (.gal) and outputs (.gif, .gal)
└── CMakeLists.txt  # Unified build configuration
```

---

## Quick Start

```bash
# 1. Build all versions
mkdir build && cd build
cmake .. && make -j

# 2. Generate initial collision data (5000 particles)
python3 scripts/generate_data.py 2500 data/inputs/g1.gal disk 0.3 0.3 0.2 0.1 0.1
python3 scripts/generate_data.py 2500 data/inputs/g2.gal disk -0.3 -0.3 0.2 -0.1 -0.1
cat data/inputs/g1.gal data/inputs/g2.gal > data/inputs/collision.gal

# 3. Run optimized simulation
./build/v5_parallel 5000 data/inputs/collision.gal 1000 0.005 8 0.5 10
mv movie.gal data/outputs/

# 4. Export visualization
python3 scripts/visualize.py data/outputs/movie.gal 5000 --save data/outputs/collision.gif --speedup 20
```

---

## Engineering Takeaways

1.  **System Bottlenecks**: Profiling revealed that malloc was the primary bottleneck, shifting focus from physics math to systems programming.
2.  **Data-Oriented Design**: Cache-friendly layout and Z-order curves demonstrate that data locality is as critical as algorithmic complexity.
3.  **Numerical Trade-offs**: Barnes-Hut involves an error margin (controlled by theta). Quantifying this error is essential for scientific software.
4.  **Robustness**: Handling edge cases like Arena Overflow and QuadTree depth limits is necessary for memory-managed systems.

---

**Deep Dives**: 
- [Technical Implementation Details](docs/TECHNICAL.md)
- [Full Performance Log](docs/BLOG.md)
