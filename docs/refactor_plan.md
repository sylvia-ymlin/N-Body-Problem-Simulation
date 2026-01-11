# Refactoring Plan: N-Body Simulation

## 1. Code Audit & Issues Analysis

### A. Memory Management (Critical)
*   **The Issue**: The current Barnes-Hut implementation uses `malloc` for every single `TNode` created during tree construction (`barnes_hut.c`: line 79) and recursively calls `free` during destruction.
*   **Impact**:
    *   **Performance**: System calls (`malloc/free`) are expensive. For $N=10,000$ particles, creating a tree might involve ~20,000 allocations *per time step*.
    *   **Fragmentation**: Severe heap fragmentation reduces cache locality.
*   **Fix**: Implement a **Linear Memory Arena (Memory Pool)**. Allocate one large block per step (or reuse one buffer) and dispense nodes linearly.

### B. Stack Safety (Critical)
*   **The Issue**: The use of Variable Length Arrays (VLA) in `galsim.c`:
    ```c
    int (*clusters)[N] = malloc(k * sizeof(*clusters));
    ```
    While `malloc` puts the data on the heap, the *type definition* involves `N` (a runtime variable), which is a C99 feature often discouraged in HPC for portability and potential compiler quirks.
    In `kmeans.c`, `int labels[N]` assigns a large array on the stack, which will cause a **Stack Overflow** if $N$ is large (e.g., > 100k, exceeding standard 8MB stack size).
*   **Fix**: Replace all VLAs with standard heap allocations (`malloc`/`calloc`) and manual index arithmetic (`idx = i * N + j`).

### C. Mathematical Efficiency
*   **The Issue**: Frequent use of `pow(dx, 2)` in Euclidean distance calculations (`kmeans.c`: line 109).
*   **Impact**: `pow()` is a general-purpose function handling logs and exponents, significantly slower than simple multiplication.
*   **Fix**: Replace with direct multiplication `dx * dx`.

### D. Build System
*   **The Issue**: The project structure has been reorganized into `src/`, `scripts/`, etc., but header includes might still rely on old paths.
*   **Fix**: Ensure `CMakeLists.txt` correctly handles include paths.

---

## 2. Refactoring Roadmap

### Phase 1: Memory Arena (High Impact)
**Goal**: Speed up Tree Construction by ~10x and eliminate recursive `free`.
1.  Define `Arena` struct in `ds.h`.
2.  Implement `arena_init`, `arena_alloc`, `arena_reset` in `barnes_hut.c`.
3.  Inject `Arena*` into `barnes_hut`, `build_tree`, and `insert` functions.
4.  Replace `malloc(sizeof(TNode))` with `arena_alloc`.
5.  Replace recursive `destroy()` with `arena_reset()`.

### Phase 2: Safety & Portability
**Goal**: Make simulation crash-proof for large $N$.
1.  Refactor `kmeans.c`: Replace `int labels[N]` with `int* labels = malloc(...)`.
2.  Refactor `galsim.c`: Replace VLA types for clusters with flat 1D arrays and manual indexing.

### Phase 3: Micro-Optimizations
**Goal**: Squeeze CPU cycles.
1.  Replace `pow(x,2)` with `x*x`.
2.  Review OpenMP scheduling (switch `static, 6` to `runtime` or `auto`).
3.  Add `const` and `restrict` keywords to pointers where applicable to help compiler vectorization.

---

## 3. Execution Plan
We will start with **Phase 1 (Memory Arena)** as it provides the most significant architectural improvement with the highest "hireability" value (shows system programming skills).
