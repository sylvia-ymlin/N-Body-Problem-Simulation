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

1. **[v1] Brute-Force Baseline**: $O(N^2)$ ground truth. Directly computed all pairwise forces with no approximation. Serves as the reference for correctness.
2. **[v2] Barnes-Hut Algorithm**: Reduced computational complexity from $O(N^2)$ to $O(N \log N)$ by approximating distant particle clusters as single nodes. This enabled scaling from thousands to tens of thousands of particles.
3. **[v3] Arena Memory Architecture**: Implemented a Linear Bump Allocator to explore scoped memory management. While experimental, it provides consistent speedups on Mac M3 architectures.
4. **[v4] Spatial Locality (Morton Curves)**: Reordered particles along a Z-order curve to improve spatial locality, reducing CPU cache misses and increasing throughput. Now supports optional Arena allocator.
5. **[v5] Parallel Load Balancing**: Leveraged OpenMP with K-Means Load Balancing to ensure uniform thread occupancy, improving parallel efficiency. Now supports optional Arena allocator.

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
