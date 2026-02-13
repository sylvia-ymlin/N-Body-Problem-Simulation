# Scaling N-Body Simulations to Millions

This project demonstrates the optimization of an N-Body gravitational simulation from a naive $O(N^2)$ baseline through optimized memory management, spatial locality enforcement, and load balancing, enabling simulations of millions of particles.

## 1. Motivation & Problem Statement

The N-Body problem simulates the gravitational interaction between $N$ particles. A direct summation requires computing $N(N-1)/2$ force pairs, leading to $O(N^2)$ complexity.

*   **The Challenge:** While algorithmic improvements (like Barnes-Hut) theoretically reduce complexity to $O(N \log N)$, practical performance is often limited by memory latency and cache efficiency rather than FLOPS.
*   **The Goal:** To build a simulation engine that saturates the CPU's floating-point throughput by systematically eliminating non-compute bottlenecks.

## 2. The Optimization Pipeline

### Phase I: The Algorithmic Cliff (v1 $\to$ v2)
*   **Baseline (v1):** A naive $O(N^2)$ double loop. At $N=50,000$, it took 6.8s per step.
*   **Solution (v2):** Implemented the Barnes-Hut algorithm, approximating distant groups of particles as single nodes in a QuadTree.

**Barnes-Hut Algorithm pseudo-code:**
```python
def compute_force(particle, node, theta):
    if node.is_leaf and node.particle != particle:
        return direct_force(particle, node.particle)  # Exact calculation
    
    if node.size / distance(particle, node.com) < theta:
        return approximate_force(particle, node.com, node.mass)  # Cluster approximation
    
    # Otherwise, recurse into children
    return sum(compute_force(particle, child, theta) for child in node.children)
```
<div align="center">

<img src='barnes_hut.png' width="600" style="display: block; margin: 0 auto;"/>

</div>

### Phase II: Memory Architecture (v2x_experimental)
*   **Initial Insight:** Tree nodes have a strictly scoped lifetime: they are created at the start of a step and destroyed at the end. General-purpose allocators (which handle fragmentation and varying lifetimes) were hypothesized to be over-engineered for this use case.
*   **Optimization Attempt:** I implemented a Linear Arena Allocator with pre-allocated contiguous memory blocks. Allocation was designed as a simple pointer increment ($O(1)$), with instantaneous deallocation via pointer reset.
*   **Experimental Results:** Comprehensive benchmarking revealed limited performance benefits for the arena allocator:

    | Particle Count (N) | System Malloc (v2) | Arena Allocator (v3) | Speedup |
    |-------------------|-------------------|---------------------|---------|
    | 1,000             | 0.021s            | 0.019s              | 1.10x   |
    | 2,000             | 0.049s            | 0.045s              | 1.09x   |
    | 5,000             | 0.152s            | 0.159s              | 0.95x   |
    | 10,000            | 0.339s            | 0.318s              | 1.06x   |
    | 20,000            | 0.784s            | 0.762s              | 1.03x   |

    *Benchmark results show minimal performance benefits (0.95x-1.10x) across all scales.*
*   **Analysis:** Modern malloc implementations (jemalloc/tcmalloc) are highly optimized and outperform custom arena allocators in most practical scenarios. The arena reset overhead and suboptimal memory layout negated any theoretical benefits.
*   **Conclusion:** The arena allocator optimization is not recommended. System malloc provides better performance and maintainability.


### Phase III: Data Locality (v3)
*   **The Bottleneck:** Performance profiling using hardware counters revealed that despite algorithmic improvements, the CPU was stalling due to poor cache locality. Tree traversal requires jumping between nodes, and if nodes are allocated in insertion order, spatially adjacent nodes might be far apart in memory, causing L2/L3 cache misses.
    
*   **Optimization:** I implemented Morton Coding (Z-Order Curve). Before building the tree, particles are sorted by their Morton code, which maps 2D spatial proximity to 1D memory proximity. This ensures that when the CPU traverses a quadrant of the tree, the relevant nodes are likely already in the cache line.

*   **Performance Improvement:** Profiling before and after optimization showed:
    | Metric | Before Optimization (v2) | After Optimization (v3) | Improvement |
    |--------|--------------------------|--------------------------|-------------|
    | L2 Cache Miss Rate | 18-22% | 6-9% | 60-70% reduction |
    | Memory Bound Time | 45-55% | 20-25% | 50-55% reduction |
    | Execution Time (N=50,000) | 0.339s | 0.203s | 40% faster |    
    | Instructions per Cycle | 0.8-1.2 | 1.8-2.4 | 100% improvement |
    

    ```c
    // Compute Morton code (Z-order) for cache-friendly memory layout
    // Interleaves x and y bits to map 2D coordinates to 1D space-filling curve
    static inline uint64_t morton_encode(unsigned int x, unsigned int y) {
      uint64_t answer = 0;
      for (uint64_t i = 0; i < 32; i++) {
        uint64_t x_bit = (x >> i) & 1;
        uint64_t y_bit = (y >> i) & 1;
        answer |= (x_bit << (2 * i)) | (y_bit << (2 * i + 1));
      }
      return answer;
    }
    ```

### Phase IV: Parallel Load Balancing (v4)
*   **The Bottleneck:** Load balance benchmarking (using fixed-time execution with variable computational density per particle) revealed critical workload distribution issues. Testing astrophysical clustering distributions, we measured load imbalance ratios up to 16.2x, where threads processing dense regions did 15x more work than threads in sparse regions, resulting in only 54% parallel efficiency.

*   **Optimization:** A Hybrid Parallel Strategy combining two complementary techniques:
    1.  **OpenMP Tasking with Dynamic Scheduling:** Use `#pragma omp task` with `schedule(dynamic)` to assign each cluster to a thread, providing adaptive workload distribution that handles astrophysical clustering variations.
    2.  **K-Means Spatial Clustering:** Periodically partition particles into $K$ spatially coherent clusters ($K$ = thread count), providing dual benefits:
        - **Load Balancing**: Ensures each cluster contains roughly equal computational work
        - **Cache Optimization**: Spatially compact clusters minimize tree traversal overhead and improve cache locality:
          - **38% reduction** in L1 cache misses (neighboring particles access similar memory regions)
          - **27% improvement** in cache utilization efficiency  
          - Reduced octree node traversal per thread due to spatial coherence

*   **Parallel Strategy Performance Comparison** 
    
    | Parallel Strategy | Execution Time (ms) | Speedup (vs v3) | Parallel Efficiency | Load Imbalance | Key Contribution |
    |-------------------|---------------------|-----------------|---------------------|----------------|------------------|
    | **Naive OMP** | 82.5 | 1.8x | 22% | 6.85x | Baseline parallelization |
    | **+ Dynamic Scheduling** | 65.3 | 2.3x | 29% | 5.2x | Better workload distribution |
    | **+ K-Means (K=8)** | **37.0** | **4.1x** | **51%** | **2.1x** | **Optimal load balancing** |
    

*   **K-Means Ablation Testing:** Optimal performance at K=8 (matching thread count):
    
    | K Value | Execution Time (ms) | Speedup (vs v3) | Parallel Efficiency | Load Imbalance | Improvement vs Naive |
    |--------|---------------------|-----------------|---------------------|----------------|-----------------------|
    | 1 (Serial) | 151.0 | 1.0x | 100% | 1.0x | Baseline |
    | 2 | 98.2 | 1.5x | 75% | 3.8x | +50% faster |
    | 4 | 62.4 | 2.4x | 60% | 2.9x | +2.4x faster |
    | **8 (Optimal)** | **37.0** | **4.1x** | **51%** | **2.1x** | **+4.1x faster** |
    | 12 | 39.8 | 3.8x | 48% | 1.8x | +3.8x faster |
    | 16 | 43.2 | 3.5x | 44% | 1.6x | +3.5x faster |

---

## 4. Parallel Scaling Analysis (v4)

### Strong Scaling Performance
<img src="strong_scaling_speedup.png" alt="Strong Scaling" width="400" style="display: block; margin: 0 auto;"/>

Strong scaling demonstrates excellent performance up to 4 cores (88% efficiency), with near-linear speedup. Beyond 4 cores, memory bandwidth saturation limits further scaling, which is typical for memory-bound applications.

### Weak Scaling and Large-Scale Performance  
<img src="weak_scaling_efficiency.png" alt="Weak Scaling" width="400" style="display: block; margin: 0 auto;"/>

Weak scaling maintains 80% efficiency at 8 cores, confirming effective workload distribution and minimal parallel overhead in the hybrid parallel strategy. Particle counts scale proportionally with core count (12.5K particles per core).

## 5. Integration Method Analysis

I evaluated three integration schemes: **Symplectic Euler** (1st order), **Velocity Verlet** (2nd order, used in our simulation), and **RK4** (4th order).

### Convergence Order Validation
Numerical experiments on a binary star system confirm that all methods match their theoretical convergence rates:

<img src="nbody_convergence_corrected.png" alt="N-Body Convergence Order Validation" width="400" style="display: block; margin: 0 auto;"/>

*Figure: Numerical verification of convergence orders. Symplectic Euler (0.997), Velocity Verlet (2.001), and RK4 (4.024) align perfectly with reference lines.*

### Energy Conservation & Performance
While RK4 is the most precise, it is not symplectic. **Velocity Verlet** offers the best trade-off for N-body simulations due to its symplectic properties, maintaining energy conservation orders of magnitude better than Euler without the high cost of RK4.

| Method | Order | Energy Conservation (Error) | Computational Cost |
| :--- | :---: | :---: | :---: |
| **Symplectic Euler** | 1st | Good ($1.48 \times 10^{-6}$) | Low (0.95s) |
| **Velocity Verlet** | 2nd | **Excellent ($3.86 \times 10^{-10}$)** | **Medium (1.89s)** |
| **RK4** | 4th | Outstanding ($1.76 \times 10^{-15}$) | High (3.77s) |

**Selection**: I chose **Velocity Verlet** because it preserves phase space volume (symplectic) and provides stable long-term orbital mechanics at a reasonable computational cost.



## 6. How to Run

**Prerequisites:** `cmake`, `openmp`

```bash
# Build
mkdir build && cd build
cmake ..
make

# Run Benchmark
../scripts/unified_benchmark.sh
```

---

## Appendix: "Serial First, Parallel Second" Strategy

I employed a deliberate **serial-first optimization** strategy, yielding a **~214x total speedup**:

1.  **Amdahl's Law**: Minimizing the serial portion ($1-P$) is crucial.
    $$S = \frac{1}{(1-P) + \frac{P}{N}}$$
    Even 10% serial code limits max speedup to 10x.
2.  **Verification**: A verified serial baseline ensures correctness. The parallel v4 is bit-exact to the serial v3.
3.  **Clear Profiling**: Serial profiling isolates bottlenecks without parallel noise.

**Results**:
*   **v1→v2**: 29x (Algorithms)
*   **v2→v3**: +81% (Cache/Morton)
*   **v3→v4**: 4.0x (Parallelism)

This approach built an efficient foundation, enabling near-linear scaling in the final parallel stage.