# Scaling N-Body Simulations to Millions

This project demonstrates the optimization of an N-Body gravitational simulation from a naive $O(N^2)$ baseline through optimized memory management, spatial locality enforcement, and load balancing, enabling simulations of millions of particles (**verified: 1M particles at 0.5 FPS on Apple M3**).

## 1. Motivation & Problem Statement

The N-Body problem simulates the gravitational interaction between $N$ particles. As particle count ($N$) grows, the computational cost explodes, making large-scale astrophysical simulations (millions of particles) computationally infeasible with brute-force methods. The primary challenge is not just reducing operation count, but overcoming the memory bandwidth bottleneck that emerges when scaling up.

*   **The Challenge:** Algorithmic improvements (like Barnes-Hut) theoretically reduce complexity to $O(N \log N)$. However, performance profiling using the Roofline Model revealed that the simulation is memory-bound, not compute-bound.

<div align="center">
  <img src="roofline_analysis.png" width="600" />
</div>

* **The Goal**: To build a simulation engine that pushes toward the hardware's limits by systematically eliminating memory latency, optimizing cache locality, and balancing parallel workloads.

<div align="center">

| Version | Breakthrough             | Same-N Speedup (N=50K) |
| :------ | :----------------------- | :--------------------- |
| **v1**  | Brute-Force $O(N^2)$     | 1x (670.31s)           |
| **v2**  | Barnes-Hut $O(N \log N)$ | ~50x (13.46s)          |
| **v3**  | Arena Allocator          | ~58x (11.58s)          |
| **v4**  | Morton + Optional Arena  | ~69x (9.73s)           |
| **v5**  | Parallel + Work-Stealing | ~149x (4.49s)          |

</div>

## 2. The Optimization Pipeline

### Phase I: Algorithmic Complexity Reduction (v1 $\to$ v2)
*   **Baseline (v1):** A naive $O(N^2)$ double loop. For $N=50,000$, it must compute $N(N-1)/2 \approx 1.25$ billion interactions per step, taking **670.31s** for just 100 steps.
*   **Solution (v2):** Implemented the Barnes-Hut algorithm, approximating distant groups of particles as single nodes in a QuadTree.
*   **Performance:** Reduced computation time from **670.31s** to **13.46s** (50x speedup) at $N=50,000$.

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

### Phase II: Memory Architecture (v3)
*   **Optimization:** A Linear Arena Allocator replaces general-purpose `malloc`. Tree nodes are allocated sequentially in contiguous memory blocks via a simple bump pointer, with $O(1)$ allocation/deallocation cost.
*   **Results:** This reduces memory overhead and improves cache locality, yielding a ~16% runtime reduction.

    | Particle Count (N) | System Malloc (v2) | Arena Allocator (v3) | Speedup |
    | ------------------ | ------------------ | -------------------- | ------- |
    | 50,000             | 13.46s             | 11.58s               | 1.16x   |

*   **Analysis:** General-purpose allocators (malloc) cause significant address discontinuity due to memory fragmentation and interspersed **metadata**. The Arena linear allocator ensures tree nodes remain contiguous during recursive insertion, increasing the probability that sibling nodes share the same **Cache Line**, thereby maximizing **Hardware Prefetcher** efficiency.

### Phase III: Data Locality (v4)
*   **Optimization:** Morton coding maps 2D spatial proximity to a 1D memory layout.
*   **Impact:** This reordering drastically reduces cross-page accesses during tree traversal, boosting **Operational Intensity**, and shifting the CPU from a memory **Latency-bound** state to a **Throughput-optimized** state.

    | Version                 | Time   | Speedup vs v2 |
    | :---------------------- | :----- | :------------ |
    | **v2 (Baseline)**       | 13.46s | 1.0x          |
    | **v4 (Morton only)**    | 10.47s | **1.29x**     |
    | **v4 (Morton + Arena)** | 9.73s  | **1.38x**     |

    *Note: Morton ordering alone provides ~29% speedup due to efficient prefetching. Combined with the Arena allocator, the total speedup reaches ~38%.*
    

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

### Phase IV: Parallel Load Balancing & Work-Stealing (v5)

**The Challenge:** In astrophysical clustering distributions, particle density differences caused a **16.2x load imbalance**, limiting parallel efficiency to 54%. We initially hypothesized that explicit **K-Means clustering** would be necessary to balance this load across threads. 

#### Ablation Study: Is K-Means Really Necessary?

To rigorously evaluate the contribution of the K-Means clustering, we conducted an ablation study isolating the parallel scheduling from the spatial clustering. 
Since `v4` already enforces spatial sorting (Morton Code), we tested whether applying OpenMP dynamic scheduling (Work-Stealing) directly to the Morton-sorted array (`k=0` mode) is sufficient, compared to the overhead of explicit K-Means clustering.

**Experimental Setup:**
- **Data**: 20,000 particles, Clustered Distribution (Disk).
- **Hardware**: 8-Core CPU (Apple M-Series).
- **Threads**: 8 threads for parallel versions.

**Results (N=20k, T=8):**

<div align="center">

| Method            | Configuration  | Parallel Strategy  | Time (10 steps) | Note                        |
| :---------------- | :------------- | :----------------- | :-------------- | :-------------------------- |
| **v4 (Baseline)** | Morton + Arena | Serial             | 0.86s           | Reference                   |
| **v5 (Ablated)**  | **k=0**        | **OpenMP Dynamic** | **0.29s**       | **Optimal (Fine-grained)**  |
| v5 (Hybrid 1:1)   | k=8            | 1 Cluster/Thread   | 1.18s           | Load Imbalance (Stragglers) |
| v5 (Hybrid 8:1)   | k=64           | 8 Clusters/Thread  | 0.48s           | Better, but has overhead    |
| v5 (Proposed)     | k=1024         | Over-decomposition | 2.13s           | Thrashing                   |

</div>

**Analysis:**
The ablation study delivers a counter-intuitive but critical insight: **Application of sophisticated clustering (K-Means) at this scale (N=20k) is counter-productive.**

1.  **k=8 (Naive 1:1 Mapping)**: Performance degrades (1.18s). With only 1 task per thread, the entire simulation waits for the thread processing the densest cluster (the "straggler problem").
2.  **k=64 (Over-decomposition)**: Performance recovers significantly (0.48s). By creating more tasks than threads (8:1 ratio), OpenMP's dynamic scheduler can balance the load—when a thread finishes a light cluster, it steals another.
3.  **k=0 (Implicit Decomposition)**: The winner (0.29s). Treating the Morton-sorted array as a source of fine-grained work allows the runtime to achieve near-perfect load balance with **zero** K-Means overhead. The spatial locality provided by Morton order is sufficient; explicit geometric clustering is redundant.

**Conclusion:**
For optimized N-Body simulations on shared-memory systems with moderate N, **Morton Coding + Dynamic Scheduling** is the "Minimal Optimal Architecture." Explicit clustering (K-Means) introduces overhead that only becomes justifiable when the cost of a "task" is massive (e.g., distributed memory transfers), not for lightweight shared-memory force calculations.

---

*   **Parallel Scaling Performance (N=50k)**
    
<div align="center">

| Threads | Time (Arena=0) | Time (Arena=1) | Speedup (vs T=1) |
| ------- | -------------- | -------------- | ---------------- |
| 1       | 8.90s          | 8.10s          | 1.0x             |
| 16      | 5.21s          | **4.49s**      | **1.8x**         |

</div>
    
    *Note: Scaling behavior at N=50k is limited by the relatively small problem size per core (~3k particles/core). However, the final speedup reaches ~149x over the naive baseline (v1).*

---

## 3. Parallel Scaling Analysis (v5)

### Strong Scaling Performance
<img src="strong_scaling_speedup.png" alt="Strong Scaling" width="400" style="display: block; margin: 0 auto;"/>

Strong scaling achieves 88% efficiency up to 4 cores, with near-linear speedup. Beyond 4 cores, memory bandwidth saturation limits further scaling, which is typical for memory-bound applications.

### Weak Scaling and Large-Scale Performance  
<img src="weak_scaling_efficiency.png" alt="Weak Scaling" width="400" style="display: block; margin: 0 auto;"/>

Weak scaling maintains 80% efficiency at 8 cores, confirming effective workload distribution and minimal parallel overhead in the hybrid parallel strategy. Particle counts scale proportionally with core count (12.5K particles per core).



## 4. How to Run

**Prerequisites:** `cmake`, `openmp`

```bash
# Build
mkdir build && cd build
cmake ..
make

# Run Benchmark
../scripts/mac_reproduce_all.sh

# Run Unified Benchmark (Linux/General)
# ../scripts/unified_benchmark.sh
```

---

## Appendix A: "Serial First, Parallel Second" Strategy

Following Amdahl's Law, the strategy prioritized single-thread algorithmic improvements (v1 $\to$ v4) to minimize the serial fraction ($1-P$). This compression of serial execution time maximizes the theoretical speedup limit for the subsequent parallel implementation (v5).

1.  **Amdahl's Law**: Minimizing the serial portion ($1-P$) is crucial.
    $$S = \frac{1}{(1-P) + \frac{P}{N}}$$
    Even 10% serial code limits max speedup to 10x.
2.  **Verification**: A verified serial baseline ensures correctness. The parallel v5 is numerically consistent with the serial v4.
3.  **Clear Profiling**: Serial profiling isolates bottlenecks without parallel noise.

**Results**:
*   **v1→v2**: ~50x (Algorithms)
*   **v2→v4**: ~1.38x (Morton + Arena)
*   **v4→v5**: ~2.16x (Parallelism)

---

## Appendix B: Integration Method Analysis

I evaluated three integration schemes: **Symplectic Euler** (1st order), **Velocity Verlet** (2nd order, used in our simulation), and **RK4** (4th order).

### Convergence Order Validation
Numerical experiments on a binary star system confirm that all methods match their theoretical convergence rates:

<img src="nbody_convergence_corrected.png" alt="N-Body Convergence Order Validation" width="400" style="display: block; margin: 0 auto;"/>

*Figure: Numerical verification of convergence orders. Symplectic Euler (0.997), Velocity Verlet (2.001), and RK4 (4.024) align perfectly with reference lines.*

### Energy Conservation & Performance
While RK4 is the most precise, it is not symplectic. **Velocity Verlet** offers the best trade-off for N-body simulations due to its symplectic properties, maintaining energy conservation orders of magnitude better than Euler without the high cost of RK4.

<div align="center">

| Method               | Order | Energy Conservation (Error) | Computational Cost |
| :------------------- | :---: | :-------------------------: | :----------------: |
| **Symplectic Euler** |  1st  |    $1.48 \times 10^{-6}$    |    Low (0.95s)     |
| **Velocity Verlet**  |  2nd  | **$3.86 \times 10^{-10}$**  | **Medium (1.89s)** |
| **RK4**              |  4th  |   $1.76 \times 10^{-15}$    |    High (3.77s)    |

</div>

**Memory Bandwidth Criticality**:
RK4 requires **4 force evaluations per time step**, which would quadruple the memory bandwidth pressure compared to a single-evaluation method. Given our Roofline analysis confirms the simulation is **Memory Bound**, increasing the memory traffic by 4x is unacceptable. **Velocity Verlet** (1 force evaluation per step) provides the optimal trade-off between symplectic accuracy and memory bandwidth efficiency.

### Simulation Parameters & Environment

**Hardware Environment**: Apple M3
**Simulation Settings**:
*   **Time Step**: $dt = 0.001$ (Velocity Verlet stability)
*   **Total Time**: $T = 0.1s$ ($100$ steps)
*   **Softening Length**: $\epsilon = 10^{-3}$ (prevents singularities)
*   **Gravitational Constant**: $G = 100/N$ (system size scaling)
