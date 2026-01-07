# Interview Prep Guide: High-Performance N-Body Simulation

## 1. Elevator Pitch (The "Hook")
"For my HPC project, I built a high-performance **N-Body Gravitational Simulation** from scratch in C. My goal was to beat the $O(N^2)$ bottleneck of naïve pairwise calculation.

I implemented the **Barnes-Hut algorithm**, which uses a Quadtree to reduce complexity to $O(N \log N)$. What made my implementation special was the use of **K-Means Clustering** to physically reorder particles in memory. This drastically reduced cache misses during tree traversal, demonstrating my ability to optimize code at the hardware architecture level."

---

## 2. Technical Deep Dive (The "Meat")

### Architecture & Optimization
*   **Barnes-Hut Algo**:
    *   *Concept*: Recursively divides space into 4 quadrants (Quadtree). Computes force between particle and tree node (center of mass) if they are sufficiently far apart ($d > L / \theta$).
    *   *Result*: Scaled to $100,000$ particles where naive methods choked at $5,000$.

*   **Memory Optimization (Cache Locality)**:
    *   *Problem*: Tree traversal jumps randomly in memory if particles are stored in random order (`pos_x[0]` is far from `pos_x[1]` in space).
    *   *Solution*: **K-Means Clustering**. I grouped physically close particles together in the array.
    *   *Impact*: Because Barnes-Hut traverses local nodes first, having neighbor particles adjacent in RAM means the CPU prefetcher loads relevant data *before* it's requested.

*   **Integrator Stability**:
    *   *Choice*: **Velocity Verlet** (Symplectic) vs Euler.
    *   *Reason*: Euler adds energy to the system (orbits spiral out). Verlet preserves the symplectic form, ensuring energy conservation for long-duration simulations.

---

## 3. "Star" Stories (Behavioral Q&A)

### Challenge: "Tell me about a bug you fought."
*   **Situation**: During parallelization with OpenMP, I noticed the simulation crashed or produced jagged orbits randomly.
*   **Root Cause**: Race conditions when multiple threads tried to update forces for the same particles, or standard `malloc` in parallel regions causing heap contention.
*   **Action**: I carefully privatized the `force` arrays and used OpenMP's `dynamic` scheduling to balance the load, as some branches of the Barnes-Hut tree are much deeper than others.
*   **Result**: Achieved correct parallel scaling without locking overhead.

### Optimization: "Why did you choose SoA over AoS?"
*   **Topic**: **Structure of Arrays (SoA)** vs **Array of Structures (AoS)**.
*   **Answer**: "I used SoA (`double pos_x[N]`) instead of AoS (`struct Particle {x, y}`). When building the tree, I only need positions. If I used AoS, loading `p.x` would also bring `p.v` into the cache line, wasting bandwidth. SoA improved bandwidth utilization during the tree construction phase."

---

## 4. Code Walkthrough Scenarios

1.  **"Show me the Barnes-Hut logic."**
    *   Open `project/barnes_hut.c`.
    *   *Point out*: The `compute_force` recursive function and the `theta` cutoff check (`width / distance < THETA_MAX`).

2.  **"Where is the Parallelism?"**
    *   Open `project/galsim.c`.
    *   *Point out*: `#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)` around the force calculation loop.

3.  **"Explain the K-Means implementation."**
    *   Open `project/kmeans.c`.
    *   *Point out*: How it iteratively updates centroids to group particles, which are then used to re-index the main arrays.

---

## 5. Key Metrics (Back of Napkin)
*   **Complexity**: $N \log N$ (verified via plotting).
*   **Speedup**: ~3.5x on a 4-core laptop (OpenMP scaling).
*   **Scale**: Capable of simulating $10^5$ particles in real-time frame rates.
