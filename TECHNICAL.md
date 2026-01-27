# Technical Deep Dive: N-Body Simulation Engine

## Architecture Overview

This project implements a high-performance gravitational N-body simulation using the Barnes-Hut algorithm with aggressive memory and cache optimizations.

### Core Components

```
src/
├── galsim.c         # Main simulation loop + integration
├── barnes_hut.c/h   # QuadTree construction & force calculation
├── morton.c/h       # Z-order curve encoding (bit-interleaving)
├── kmeans.c/h       # Spatial clustering for load balancing
└── ds.h             # Linear arena allocator + data structures

scripts/
├── bench_alloc.c        # Memory allocator benchmark
├── code_complexity.py   # Complexity verification plots
└── generate_data.py     # Initial condition generator
```

## Algorithm Details

### Barnes-Hut Treecode

Hierarchical space partitioning using a **QuadTree** (2D) to approximate distant particle clusters as single point masses.

**Theta Criterion**: If $\frac{s}{d} < \theta$ (where $s$ = cell width, $d$ = distance), treat entire cell as one mass.

**Complexity Analysis**:
- **Tree Construction**: $O(N \log N)$
  - Inserting one particle into a QuadTree of depth $d$ takes $O(d)$
  - For uniform distribution, $d \approx \log_4 N$
  - Total: $N \times O(\log N) = O(N \log N)$
  
- **Force Calculation**: $O(N \log N)$
  - Each particle traverses tree to depth $\approx \log N$
  - Theta criterion prunes branches, keeping traversal logarithmic
  
- **Total per timestep**: $O(N \log N)$ vs $O(N^2)$ for brute force

**Accuracy**: Configurable via $\theta$ parameter
- $\theta = 0.0$: Exact $O(N^2)$ (no approximation)
- $\theta = 0.5$: Balanced (default, <2% error)
- $\theta = 1.0$: Fast but less accurate

### Spatial Data Structures

#### Z-Order Curve (Morton Codes)

Maps 2D coordinates to 1D index while preserving spatial locality.

**Implementation**: 64-bit bit-interleaving
```c
// src/morton.c
static inline uint64_t split_by_3(unsigned int a) {
    uint64_t x = a & 0x1fffff;  // Mask to 21 bits
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8)  & 0x100f00f00f00f00f;
    x = (x | x << 4)  & 0x10c30c30c30c30c3;
    x = (x | x << 2)  & 0x1249249249249249;
    return x;
}

uint64_t morton_encode(unsigned int x, unsigned int y) {
    return split_by_3(x) | (split_by_3(y) << 1);
}
```

**Why it works**: Particles close in 2D space get close indices in 1D array → sequential memory access → better cache utilization.

#### K-Means Clustering

Pre-partitions particles into spatial clusters for OpenMP thread distribution.

**Purpose**:
- Ensures balanced workloads (each cluster has similar particle count)
- Improves cache locality (threads work on spatially coherent regions)
- Acts as pre-conditioner for Z-order sorting

**Algorithm**: Lloyd's algorithm with spatial centroids
- Initialize $k$ centroids (typically $k$ = number of threads)
- Iterate: assign particles to nearest centroid, recompute centroids
- Convergence: typically 5-10 iterations

## Memory Management

### Arena Allocator Implementation

**Problem**: Standard `malloc/free` for tree nodes causes:
- Syscall overhead (kernel context switches)
- Heap fragmentation (non-contiguous allocations)
- Poor cache locality (random memory layout)

**Solution**: Linear arena allocator (bump pointer)

```c
// src/ds.h
typedef struct NodeArena {
    TNode *buffer;    // Pre-allocated contiguous block
    size_t capacity;  // Total nodes available
    size_t used;      // Current allocation pointer
} NodeArena;

// O(1) allocation - just increment pointer
TNode* arena_alloc(NodeArena* arena) {
    if (arena->used >= arena->capacity) {
        fprintf(stderr, "Arena exhausted\n");
        exit(1);
    }
    return &arena->buffer[arena->used++];
}

// O(1) bulk deallocation - reset pointer
void arena_reset(NodeArena* arena) {
    arena->used = 0;  // All nodes "freed" instantly
}
```

**Key insight**: Tree nodes have identical lifetimes—all created at step start, all destroyed at step end. No need for individual `free()`.

**Memory layout benefits**:
- Nodes allocated sequentially in memory
- Tree traversal accesses contiguous memory regions
- CPU prefetcher can predict access patterns
- L1/L2 cache hit rates dramatically improved

### Benchmark Results

**Test setup**: 100 iterations of allocating/freeing 100,000 tree nodes

| Allocator | Time | Operations/sec | Speedup |
|-----------|------|----------------|---------|
| System `malloc`/`free` | 0.42s | 23.8M ops/s | 1x |
| **Linear Arena** | **< 0.001s** | **> 10B ops/s** | **> 400x** |

**Why such massive speedup?**
- `malloc`: Traverses free lists, updates metadata, syscalls
- Arena: Single pointer increment (1 CPU instruction)

## Numerical Integration

### Symplectic Euler (First-Order)
```c
// Update velocity, then position
vx[i] += fx[i] * dt;
vy[i] += fy[i] * dt;
px[i] += vx[i] * dt;
py[i] += vy[i] * dt;
```
- Fast but accumulates energy error
- Suitable for short simulations

### Velocity Verlet (Second-Order)
```c
// Half-step velocity update
vx[i] += 0.5 * fx[i] * dt;
vy[i] += 0.5 * fy[i] * dt;

// Full-step position update
px[i] += vx[i] * dt;
py[i] += vy[i] * dt;

// Recompute forces, then complete velocity update
compute_forces();
vx[i] += 0.5 * fx_new[i] * dt;
vy[i] += 0.5 * fy_new[i] * dt;
```
- Energy-conserving (symplectic)
- Default choice for long simulations

### Plummer Softening

Prevents singularities when particles are very close:

$$F = \frac{G \cdot m_1 \cdot m_2}{(r^2 + \epsilon^2)^{3/2}} \cdot \vec{r}$$

where $\epsilon$ is the softening parameter (typically $\epsilon \approx 0.01$).

## Parallelization Strategy

### OpenMP Implementation

```c
// Force calculation: Read-only tree access → trivially parallel
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
for (int i = 0; i < N; i++) {
    double fx_local = 0.0, fy_local = 0.0;
    barnes_hut_force(particles[i], tree_root, &fx_local, &fy_local);
    fx[i] = fx_local;
    fy[i] = fy_local;
}
```

**Key design choices**:

1. **Dynamic scheduling**: `schedule(dynamic, CHUNK_SIZE)`
   - Tree depth varies with particle density
   - Static scheduling causes thread starvation
   - Dynamic assigns work on-demand

2. **No race conditions**:
   - Tree is read-only during force calculation
   - Each thread writes to disjoint array indices
   - No locks or atomics needed

3. **K-Means pre-partitioning**:
   - Ensures threads work on spatially coherent regions
   - Maximizes L2 cache reuse within threads
   - Reduces false sharing

### Load Balancing Analysis

**Without K-Means**: Threads assigned sequential particle ranges
- Thread 1 might get dense cluster (deep tree, slow)
- Thread 2 might get sparse region (shallow tree, fast)
- Result: Thread 1 becomes bottleneck

**With K-Means**: Threads assigned spatial clusters
- Each cluster has similar particle density
- Tree depth roughly equal across clusters
- Result: Balanced workload, near-linear scaling

## Benchmarking Methodology

### Tree Traversal Optimization

Comparing recursive vs stackless iterative traversal:

**Recursive implementation**:
```c
void traverse_recursive(TNode* node, Particle p) {
    if (is_leaf(node)) {
        compute_force(p, node);
        return;
    }
    if (theta_criterion(node, p)) {
        compute_force(p, node);  // Approximate
    } else {
        for (int i = 0; i < 4; i++) {
            if (node->children[i])
                traverse_recursive(node->children[i], p);
        }
    }
}
```

**Stackless implementation**:
```c
void traverse_stackless(TNode* root, Particle p) {
    TNode* stack[MAX_DEPTH];
    int top = 0;
    stack[top++] = root;
    
    while (top > 0) {
        TNode* node = stack[--top];
        if (is_leaf(node) || theta_criterion(node, p)) {
            compute_force(p, node);
        } else {
            for (int i = 0; i < 4; i++) {
                if (node->children[i])
                    stack[top++] = node->children[i];
            }
        }
    }
}
```

**Results** ($N=100k$, 200 timesteps):
| Method | Time | Speedup |
|--------|------|---------|
| Recursive | 0.0619s | 1x |
| **Stackless** | **0.0555s** | **1.11x** |

**Why stackless is faster**:
- Eliminates function call overhead (register save/restore)
- Better instruction cache locality (single function body)
- Explicit stack control (can optimize layout)

### Accuracy Verification

Comparing against reference brute-force implementation:

<p align="center">
  <img src="docs/code_complexity.png" width="500" />
</p>

**Test setup**: 1000-particle elliptical galaxy, 200 timesteps

| Theta | Relative Error | Speedup vs Brute Force |
|-------|----------------|------------------------|
| 0.0 | 0% (exact) | 1x |
| 0.3 | <0.5% | ~15x |
| 0.5 | <2% | ~30x |
| 0.7 | <5% | ~50x |
| 1.0 | <10% | ~80x |

**Conclusion**: $\theta = 0.5$ provides excellent accuracy/speed balance.

## Engineering Decisions

### Why C over C++?

**Explicit memory control**:
- No hidden allocations (vtables, RTTI, exceptions)
- Transparent pointer arithmetic for arena allocator
- Easier to reason about cache behavior

**Performance transparency**:
- No constructor/destructor overhead
- No virtual function indirection
- Direct mapping to assembly

**Simplicity**:
- Smaller language surface area
- Easier to audit for performance issues
- Better for teaching low-level optimization

### Why QuadTree over Octree?

**2D simulation benefits**:
- Reduces tree depth (4-way vs 8-way branching)
- Simpler to visualize and debug
- Faster tree construction

**Extensibility**:
- Can extend to 3D with minimal changes
- Same algorithmic principles apply
- Good starting point for learning

### Why Linear Arena over Pool Allocator?

**Lifetime analysis**:
- All tree nodes created at timestep start
- All destroyed at timestep end
- No need for individual deallocation

**Performance**:
- Pool allocator: maintains free list, requires bookkeeping
- Linear arena: single pointer increment
- Arena is simpler and faster for this use case

**Cache locality**:
- Linear arena guarantees contiguous allocation
- Pool allocator may fragment over time
- Contiguous nodes → better prefetching

### Why SoA over AoS?

**Structure of Arrays (SoA)** - Current approach:
```c
double pos_x[N];
double pos_y[N];
double vel_x[N];
double vel_y[N];
```

**Array of Structures (AoS)** - Alternative:
```c
struct Particle {
    double pos_x, pos_y;
    double vel_x, vel_y;
} particles[N];
```

**Why SoA wins**:
- SIMD-friendly: Can process 4 x-coordinates at once
- Cache-efficient: Only load needed fields (e.g., positions for tree building)
- Better memory bandwidth utilization

## Performance Tuning Guide

### Accuracy vs Speed

Adjust the `theta` parameter:
```bash
./galsim N input.gal output.gal steps theta
```

- `theta = 0.0`: Exact $O(N^2)$ (no approximation)
- `theta = 0.5`: Balanced (default, <2% error)
- `theta = 1.0`: Fast but less accurate (~10% error)

### Parallelism

Set OpenMP threads via environment variable:
```bash
export OMP_NUM_THREADS=8
./galsim 10000 input.gal output.gal 200 0.5
```

**Optimal thread count**: Typically number of physical cores (not hyperthreads).

### Memory Tuning

For very large simulations, adjust arena size in `ds.h`:
```c
#define ARENA_SIZE (3 * N)  // Increase multiplier if tree is deep
```

**Rule of thumb**: 
- Uniform distribution: 2-3x particles
- Clustered distribution: 4-5x particles

### Compiler Optimization

Build with maximum optimization:
```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

**Key flags** (automatically set by CMake):
- `-O3`: Aggressive optimization
- `-march=native`: Use CPU-specific instructions
- `-fopenmp`: Enable OpenMP

## Scalability Analysis

### Weak Scaling

Fixed particles per thread, increase thread count:

| Threads | Particles | Time | Efficiency |
|---------|-----------|------|------------|
| 1 | 10k | 1.0s | 100% |
| 2 | 20k | 1.1s | 91% |
| 4 | 40k | 1.2s | 83% |
| 8 | 80k | 1.4s | 71% |

**Observation**: Good scaling up to 4 threads, then overhead increases.

### Strong Scaling

Fixed problem size (100k particles), increase threads:

| Threads | Time | Speedup | Efficiency |
|---------|------|---------|------------|
| 1 | 10.0s | 1.0x | 100% |
| 2 | 5.3s | 1.9x | 95% |
| 4 | 2.8s | 3.6x | 90% |
| 8 | 1.6s | 6.3x | 79% |

**Observation**: Near-linear scaling up to 4 threads, diminishing returns beyond.

## Future Optimizations

### SIMD Vectorization

Explicitly vectorize force calculations using AVX2/AVX-512:
```c
// Process 4 particles at once
__m256d px_vec = _mm256_load_pd(&pos_x[i]);
__m256d py_vec = _mm256_load_pd(&pos_y[i]);
// ... vectorized distance calculation
```

**Expected gain**: 2-4x for force calculation kernel.

### GPU Acceleration

**Challenges**:
- Tree traversal is irregular (branch divergence)
- Recursive algorithms don't map well to SIMD

**Potential solutions**:
- Flatten tree into array (BFS layout)
- Use warp-level primitives for divergence
- Hybrid: CPU builds tree, GPU computes forces

**Expected gain**: 10-50x for large $N$ (>1M particles).

### MPI for Multi-Node

For scaling beyond single machine ($N > 10^9$):
- Domain decomposition (Orthogonal Recursive Bisection)
- Halo exchange for boundary particles
- Tree replication or distributed tree construction

**Expected gain**: Linear scaling up to 100s of nodes.

---

## References

- Barnes, J., & Hut, P. (1986). "A hierarchical O(N log N) force-calculation algorithm". *Nature*, 324(6096), 446-449.
- Warren, M. S., & Salmon, J. K. (1993). "A parallel hashed Oct-Tree N-body algorithm". *Supercomputing*, 12-21.
- Dehnen, W. (2002). "A Hierarchical O(N) Force Calculation Algorithm". *Journal of Computational Physics*, 179(1), 27-42.
