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
2. **[v2] Barnes-Hut Breakthrough**: Switched to hierarchical space partitioning (Barnes-Hut, $O(N \log N)$). Runtime was limited by memory allocator overhead.
3. **[v2x] Arena Memory Architecture (Experimental)**: Replaced `malloc` with a Linear Bump Allocator ([implementation](v2x_arena_experimental/ds.h)) to explore memory management optimizations. While theoretically faster for scoped allocations, performance analysis showed limited benefits compared to modern malloc implementations in practical scenarios.
4. **[v3] Spatial Locality (Morton Curves)**: Reordered particles along a Z-order curve ([implementation](v3_morton/morton.c)) to improve spatial locality, reducing CPU cache misses and increasing throughput.
5. **[v4] Parallel Load Balancing**: Leveraged OpenMP with K-Means Load Balancing ([implementation](v4_parallel/kmeans.c)) to ensure uniform thread occupancy even in highly non-uniform star clusters, improving parallel efficiency.

| Version   | Breakthrough             | Max Scale | Same-N Speedup (N=50K) |
| :-------- | :----------------------- | :-------- | :--------------------- |
| **v1**    | Brute-Force $O(N^2)$     | 4.2K      | 1x (66.5s)             |
| **v2**    | Barnes-Hut $O(N \log N)$ | 50K       | ~29x (2.3s)            |
| **v3**    | Memory/Cache Optimized   | 500K      | ~53x (1.2s)            |
| **v4**    | Parallel Engine          | 2.5M      | ~214x (0.31s)          |

**Numerical Stability**: Supports Velocity Verlet (energy conservation) and RK4 (high-precision inquiry).



## Project Structure

```text
.
├── v1_naive: Naive $O(N^2)$ implementation. (Verified Correct)
├── v2_barnes_hut: Basic Barnes-Hut implementation. (Verified Correct)
├── v2x_arena_experimental: Experimental Arena allocator for memory management
├── v3_morton: Optimized memory layout using Z-order (Morton) sorting. (Verified Correct)
├── v4_parallel: Parallel implementation using OpenMP + K-Means Load Balancing. (Verified Correct)
├── scripts/            # Unified Benchmark & Validation Suite
├── data/               # Structured inputs/outputs
├── docs/               # Technical Report & Performance Logs
```
