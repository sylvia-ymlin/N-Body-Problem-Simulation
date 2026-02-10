# N-Body Simulation Engine: Scaling to Millions on a Consumer PC

<p align="center">
  <img src="https://img.shields.io/badge/Language-C-00599C" alt="C">
  <img src="https://img.shields.io/badge/Parallelism-OpenMP-blue" alt="OpenMP">
  <img src="https://img.shields.io/badge/Scale-2.5M_Particles-green" alt="Scale">
</p>

**"How can we leverage algorithmic precision and hardware-aware systems engineering to run million-scale astronomical simulations efficiently on a standard laptop?"**

This project documents the journey from a 4,200-particle brute-force simulation to a 2.5-million-particle parallel engine—an ~80x same-scale speedup and a 600x increase in maximum tractable problem size. It is an engineering case study in identifying system bottlenecks and architecting custom memory/cache solutions to break through them.

![Cinematic simulation of a galactic collision](data/outputs/cinematic_collision.gif)

---

## The Optimization Story

My objective was to push the simulation boundary from 4,200 particles to 2.5 million through five evolutionary stages:

1.  **[v1] Brute-Force Baseline**: Established $O(N^2)$ ground truth. Identified that physics math isn't the bottleneck—algorithmic complexity is.
2.  **[v2] Barnes-Hut Breakthrough**: Switched to hierarchical space partitioning ($O(N \log N)$). Paradoxically, performance hit a wall due to System Allocator Overhead.
3.  **[v3] Arena Memory Architecture**: Replaced `malloc` with a Linear Bump Allocator (see [implementation](v3_arena/ds.h)). Achieved a 400x reduction in allocation time by aligning memory lifetimes with simulation steps.
4.  **[v4] Spatial Locality (Morton Curves)**: Addressed CPU cache misses. By reordering particles along a Z-order curve (see [implementation](v4_morton/morton.c)), I ensured that spatial neighbors are also memory neighbors.
5.  **[v5] Parallel Load Balancing**: Leveraged OpenMP with K-Means Load Balancing (see [implementation](v5_parallel/kmeans.c)) to ensure uniform thread occupancy even in highly non-uniform star clusters.

### Key Results
| Version   | Breakthrough             | Max Scale | Same-N Speedup (N=50K) |
| :-------- | :----------------------- | :-------- | :--------------------- |
| **v1**    | Brute-Force $O(N^2)$     | 4.2K      | 1x (67.1s)             |
| **v2**    | Barnes-Hut $O(N \log N)$ | 50K       | ~22x (3.1s)            |
| **v3/v4** | Memory/Cache Optimized   | 500K      | ~26x (2.5s)            |
| **v5**    | Parallel Engine          | 2.5M      | ~83x (0.81s)           |

---

## Technical Exhibit

For a comprehensive engineering report detailing the specific optimizations, memory architecture, and validation results, see [REPORT.md](docs/REPORT.md).

- **Numerical Stability**: Supports Velocity Verlet (energy conservation) and RK4 (high-precision inquiry).
- **Verification**: 
  - `validate_accuracy.py`: Quantifies approximation errors against the $O(N^2)$ baseline.
  - `unified_benchmark.sh`: Standardized rig to reproduce optimization claims.

---

## Quick Start

Build and run the 2.5M particle simulation:

```bash
./scripts/unified_benchmark.sh
```

---

## Project Structure

```text
.
├── v1_naive: Naive $O(N^2)$ implementation. (Verified Correct)
├── v2_barnes_hut: Basic Barnes-Hut implementation. (Verified Correct)
├── v3_arena: Optimized memory management using Arena allocator. (Verified Correct)
├── v4_morton: Optimized memory layout using Z-order (Morton) sorting. (Verified Correct)
├── v5_parallel: Parallel implementation using Pthreads/OpenMP + K-Means Load Balancing. (Verified Correct)
├── scripts/            # Unified Benchmark & Validation Suite
├── data/               # Structured inputs/outputs
├── docs/               # Technical Report & Performance Logs
```
