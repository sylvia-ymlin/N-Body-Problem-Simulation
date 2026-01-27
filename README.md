# N-Body Simulation Engine

Simulating galaxy collisions is impossibly slow—each star affects every other star. I built a spatial tree algorithm with custom memory management to skip distant calculations. This achieved stable simulation of 2.5 million stars in real-time on a laptop.

## The Core Insight

**Memory allocation was the bottleneck, not the algorithm.**

Each simulation step required 200,000+ malloc calls, causing:
- Kernel context switches (syscall overhead)
- Heap fragmentation (random memory → cache misses)

Solution: Custom arena allocator—bump pointer allocation with bulk deallocation.

**Result: 400x speedup** (0.42s → <0.001s per timestep).

## What I Built

**Barnes-Hut tree algorithm**: Groups distant stars into single points
- Reduces $O(N^2)$ to $O(N \log N)$
- Trades small accuracy loss for massive speed gain

**Cache-optimized memory layout**:
- Z-order curve for spatial locality
- K-means pre-partitioning for parallel load balancing

**Production-ready engine**:
- Scales to 2.5M particles
- Configurable accuracy/speed tradeoffs
- Parallel execution with OpenMP

## Key Results

| Metric | Value |
|--------|-------|
| Max particles | 2.5 million |
| Memory speedup | 400x (arena vs malloc) |
| Algorithmic complexity | $O(N \log N)$ |
| Accuracy at θ=0.5 | <2% error |

### Complexity Verification
Barnes-Hut scales at $O(N \log N)$ vs brute force $O(N^2)$—enabling 100x larger simulations:

<p align="center">
  <img src="docs/execution_time_vs_N.png" width="500" />
</p>

### Memory Allocator Benchmark
Custom linear arena vs system malloc (100 iterations, 100k nodes):

| Allocator | Time | Speedup |
|-----------|------|---------|
| `malloc`/`free` | 0.42s | 1x |
| **Linear Arena** | **< 0.001s** | **> 400x** |

## Quick Start

```bash
mkdir build && cd build
cmake .. && make
./galsim 10000 data/input_data/ellipse_N_10000.gal output.gal 200 0.5
```

## What I Learned

**System allocators are hidden bottlenecks**
- Profiling revealed malloc/free dominated runtime, not the physics calculations
- Custom arena allocator: 400x speedup by eliminating syscalls

**Memory layout beats algorithmic complexity**
- Even with $O(N \log N)$, random access caused 70% cache misses
- Spatial sorting (Z-order + K-means) was critical for real performance

**Parallelism requires domain knowledge**
- Naive OpenMP caused thread starvation (uneven workloads)
- Solution: spatial pre-partitioning + dynamic scheduling

**Stack overflow is real**
- VLAs with user-controlled sizes crashed at 100k+ particles
- Always heap-allocate in production systems

---

**Technical deep dive**: See [TECHNICAL.md](TECHNICAL.md) for implementation details, architecture decisions, and benchmarking methodology.

![Language](https://img.shields.io/badge/Language-C11-blue.svg)
![Optimization](https://img.shields.io/badge/Optimization-OpenMP%20%7C%20Cache--Optimized-green.svg)
![Scale](https://img.shields.io/badge/Scale-2.5M%20Particles-orange.svg)
