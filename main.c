#include "simulation.h"
#include "core/io.h"
#include "core/time_utils.h"
#include "core/types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  if (argc < 9) {
    printf("N-Body Simulation Engine: Scalable Optimization Path\n");
    printf("Usage: %s <version> N <input.gal> nsteps dt n_threads theta k\n",
           argv[0]);
    printf("Versions:\n");
    printf("  1: Naive O(N²) - Reference baseline\n");
    printf("  2: Barnes-Hut   - O(N log N) algorithm\n");
    printf("  3: Arena        - Memory prefetching optimization\n");
    printf("  4: Morton       - Cache locality (Z-order) optimization\n");
    printf("  5: Parallel     - Scaling (OpenMP Dynamic + Morton)\n");
    return 1;
  }

  int version = atoi(argv[1]);
  int N = atoi(argv[2]);
  char *filename = argv[3];
  int nsteps = atoi(argv[4]);
  double dt = atof(argv[5]);
  int n_threads = atoi(argv[6]);
  double theta = atof(argv[7]);
  int k_clusters = atoi(argv[8]);

  KernelConfig config = {theta, n_threads, k_clusters};

  ParticleSystem sys = io_read_particles(filename, N);

  ForceComputeKernel kernel = NULL;
  switch (version) {
  case 1:
    kernel = compute_force_v1_naive;
    break;
  case 2:
    kernel = compute_force_v2_barnes_hut;
    break;
  case 3:
    kernel = compute_force_v3_arena;
    break;
  case 4:
    kernel = compute_force_v4_morton;
    break;
  case 5:
    kernel = compute_force_v5_parallel;
    break;
  default:
    fprintf(stderr, "Version %d not found\n", version);
    return 1;
  }

  printf("🚀 Running v%d | N=%d | Steps=%d | Threads=%d\n", version, N, nsteps,
         n_threads);
  double start = sim_time_now();

  // Initial Force
  kernel(&sys, &config);

  // Main Simulation Loop
  for (int step = 0; step < nsteps; step++) {
    // 1. Kick 1 & Drift (Velocity Verlet)
    for (int i = 0; i < N; i++) {
      double m_inv = 1.0 / sys.mass[i];
      sys.vx[i] += 0.5 * dt * sys.fx[i] * m_inv;
      sys.vy[i] += 0.5 * dt * sys.fy[i] * m_inv;
      sys.pos_x[i] += dt * sys.vx[i];
      sys.pos_y[i] += dt * sys.vy[i];
    }

    // 2. Recompute Forces
    kernel(&sys, &config);

    // 3. Kick 2
    for (int i = 0; i < N; i++) {
      double m_inv = 1.0 / sys.mass[i];
      sys.vx[i] += 0.5 * dt * sys.fx[i] * m_inv;
      sys.vy[i] += 0.5 * dt * sys.fy[i] * m_inv;
    }

    if (step % 50 == 0)
      printf("Step %d/%d\r", step, nsteps);
    fflush(stdout);
  }

  double end = sim_time_now();
  printf("\n✅ Simulation Complete: %.2fs\n", end - start);

  char out_name[64];
  sprintf(out_name, "data/outputs/result_v%d.gal", version);
  io_write_result(out_name, &sys);
  io_free_particles(&sys);
  return 0;
}
