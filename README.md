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
5. **[v5] Parallel Load Balancing**: Mitigates severe astrophysical load imbalances by decoupling spatial locality from computational weight. Uses K-Means clustering paired with OpenMP `dynamic` scheduling for efficient work-stealing.

| Version   | Breakthrough             | Same-N Speedup (N=50K) |
| :-------- | :----------------------- | :--------------------- |
| **v1**    | Brute-Force $O(N^2)$     | 1x (670.31s)           |
| **v2**    | Barnes-Hut $O(N \log N)$ | ~50x (13.46s)          |
| **v3**    | Arena Allocator          | ~58x (11.58s)          |
| **v4**    | Morton + Optional Arena  | ~69x (9.73s)           |
| **v5**    | Parallel + Optional Arena| ~149x (4.49s)          |

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
├── v1_naive: Naive $O(N^2)$ implementation. (Verified Correct)
├── v2_barnes_hut: Basic Barnes-Hut implementation. (Verified Correct)
├── v3_arena: Arena allocator for memory management. (Verified Correct)
├── v4_morton: Optimized memory layout using Z-order (Morton) sorting + Optional Arena. (Verified Correct)
├── v5_parallel: Parallel implementation using OpenMP + K-Means Load Balancing + Optional Arena. (Verified Correct)
├── scripts/            # Unified Benchmark & Validation Suite
├── data/               # Structured inputs/outputs
├── docs/               # Technical Report & Performance Logs
```
