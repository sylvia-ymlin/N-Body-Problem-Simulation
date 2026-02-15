# N-Body Simulation Engine: Scaling to Millions on a Consumer PC

<p align="center">
  <img src="https://img.shields.io/badge/Language-C-00599C" alt="C">
  <img src="https://img.shields.io/badge/Parallelism-OpenMP-blue" alt="OpenMP">
</p>

<div style="text-align: center;">
<img src="docs/cinematic_collision.gif" width="400" alt="Galactic Collision GIF">
</div>

  <br>

This project iteratively extends and accelerates the simulation engine, pushing the boundary from thousands to millions of particles. For a comprehensive engineering report detailing the specific optimizations, memory architecture, and validation results, see [REPORT.md](docs/REPORT.md).

1. **[v1] Brute-Force Baseline**: $O(N^2)$ ground truth. Directly computes all pairwise forces with no approximation. Serves as the numerical reference.
2. **[v2] Barnes-Hut Algorithm**: Reduces computational complexity to $O(N \log N)$ by approximating distant particle clusters via a QuadTree. Enables scaling to tens of thousands of particles.
3. **[v3] Arena Memory Architecture**: Replaces `malloc` with an $O(1)$ Linear Bump Allocator. Eliminates metadata fragmentation and aligns tree nodes in memory to maximize Hardware Prefetcher efficiency.
4. **[v4] Spatial Locality (Morton Curves)**: Maps 2D coordinates to a 1D Z-order curve before tree construction. This enforces strict memory contiguity, boosting Operational Intensity and shifting the CPU towards throughput-optimization.
5. **[v5] Parallel Load Balancing**: Implements fine-grained task scheduling to mitigate severe astrophysical load imbalances. **Ablation studies** reveal that while K-Means clustering was initially explored for geometric partitioning, the most efficient architecture utilizes **OpenMP dynamic scheduling** directly over **Morton-ordered arrays**, achieving near-perfect load balance with zero clustering overhead.

<div align="center">

| Version | Breakthrough                 | Same-N Speedup (N=50K) |
| ------- | ---------------------------- | ---------------------- |
| **v1**  | Brute-Force $O(N^2)$         | 1x (670.31s)           |
| **v2**  | Barnes-Hut $O(N \log N)$     | ~50x (13.46s)          |
| **v3**  | Arena Allocator              | ~58x (11.58s)          |
| **v4**  | Morton + Optional Arena      | ~69x (9.73s)           |
| **v5**  | **Parallel + Work-Stealing** | **~149x (4.49s)**      |

</div>

## Engineering Highlights & Critical Insights

### 1. The "Million Particle" Milestone
The engine has been verified to simulate **1 Million particles** at **0.5 FPS** on a consumer-grade **Apple M3** (8-core). This demonstrates the power of shifting a system from a memory-latency-bound state to a throughput-optimized state.

### 2. Lessons from the Ablation Study: Why K-Means Failed
One of the most significant findings of this project was the dismissal of explicit geometric clustering (K-Means).
* **Redundancy**: Morton coding already provides sufficient spatial locality for cache efficiency.
* **Density bottlenecks**: Geometric partitioning cannot resolve computational density spikes (e.g., galactic cores).
* **The Solution**: Leveraging the **Work-Stealing** nature of OpenMP's `dynamic` scheduler over a 1D space-filling curve (Z-order) is the "Minimal Optimal Architecture" for this problem.

### 3. Memory Bandwidth vs. Numerical Precision
While higher-order methods like RK4 offer better single-step precision, the **Roofline Model** confirms that the increased memory traffic (4 evaluations per step) would saturate the bandwidth. **Velocity Verlet** was selected as the optimal compromise between symplectic energy conservation and hardware memory limits.

## Quick Start

```bash
# 1. Build
mkdir build && cd build
cmake .. && make

# 2. Tests all versions
../scripts/mac_reproduce_all.sh
```

## Project Structure

```text
.
├── v1_naive/           # Baseline: Naive O(N²) implementation
├── v2_barnes_hut/      # Algorithmic shift: QuadTree-based O(N log N)
├── v3_arena/           # Memory optimization: Linear Arena Allocator
├── v4_morton/          # Cache optimization: Z-order (Morton) sorting
├── v5_parallel/        # Final parallel version. Supports two modes:
│                       #  - Mode 1: Morton + OpenMP Dynamic (The Winner)
│                       #  - Mode 2: K-Means Load Balancing (Experimental)
├── scripts/            # Python-based benchmarking, validation & plotting
├── data/               # Input profiles (.gal) and simulation outputs
├── docs/               # Technical Report, Roofline analysis & diagrams
├── CMakeLists.txt      # Unified build system
└── README.md           # Project overview and summary
```
