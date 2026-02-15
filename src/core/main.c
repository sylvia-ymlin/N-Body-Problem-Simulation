#include "../kernels/kernels.h"
#include "io.h"
#include "types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double get_time() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
#endif
}

int main(int argc, char *argv[]) {
  // Usage: ./nbody_simulate <version> N filename nsteps delta_t n_threads
  // theta_max k [use_arena]
  if (argc < 9) {
    printf("Usage: %s <version> N filename nsteps delta_t n_threads theta_max "
           "k [use_arena]\n",
           argv[0]);
    printf("Versions:\n");
    printf("  1: Naive O(N^2)\n");
    printf("  2: Barnes-Hut (Malloc)\n");
    printf("  3: Barnes-Hut (Arena)\n");
    printf("  4: Barnes-Hut (Morton)\n");
    printf("  5: Parallel Barnes-Hut\n");
    return 1;
  }

  int version = atoi(argv[1]);
  int N = atoi(argv[2]);
  char *filename = argv[3];
  int nsteps = atoi(argv[4]);
  double dt = atof(argv[5]);
  int n_threads = atoi(argv[6]);
  double theta_max = atof(argv[7]);
  int k_clusters = atoi(argv[8]);
  int use_arena = 0;
  if (argc > 9)
    use_arena = atoi(argv[9]);

  KernelConfig config;
  config.theta_max = theta_max;
  config.n_threads = n_threads;
  config.use_arena = use_arena;
  config.k_clusters = k_clusters;

  // Load Particles
  ParticleSystem sys = io_read_particles(filename, N);

  // Select Kernel
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
    fprintf(stderr, "Invalid version %d\n", version);
    io_free_particles(&sys);
    return 1;
  }

  // Initialize/Clear movie file
  FILE *movie = fopen("movie.gal", "wb");
  if (movie)
    fclose(movie);

  printf("Starting simulation version %d with N=%d, steps=%d, threads=%d\n",
         version, N, nsteps, n_threads);
  double start_time = get_time();

  // 1. Initial Force Compute
  kernel(&sys, &config);

  // 2. Integration Loop
  for (int step = 0; step < nsteps; step++) {
    // Output frame
    if (step % 10 == 0) {
      io_write_frame("movie.gal", &sys);
    }

    // Kick 1 & Drift
    for (int i = 0; i < N; i++) {
      double ax = sys.fx[i] / sys.mass[i];
      double ay = sys.fy[i] / sys.mass[i];
      sys.vx[i] += 0.5 * dt * ax;
      sys.vy[i] += 0.5 * dt * ay;
      sys.pos_x[i] += dt * sys.vx[i];
      sys.pos_y[i] += dt * sys.vy[i];
    }

    // Recompute Force
    kernel(&sys, &config);

    // Kick 2
    for (int i = 0; i < N; i++) {
      double ax = sys.fx[i] / sys.mass[i];
      double ay = sys.fy[i] / sys.mass[i];
      sys.vx[i] += 0.5 * dt * ax;
      sys.vy[i] += 0.5 * dt * ay;
    }
  }

  double elapsed = get_time() - start_time;
  printf("Simulation took %f seconds\n", elapsed);

  io_write_result("result.gal", &sys);
  io_free_particles(&sys);

  return 0;
}
