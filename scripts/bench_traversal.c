#include "../src/barnes_hut.h"
#include "../src/ds.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Mock implementation of helper functions needed if linking fails due to
// missing galsim dependencies But since we link against barnes_hut.c, we should
// be fine.

double get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

int main(int argc, char **argv) {
  int N = 100000;
  if (argc > 1)
    N = atoi(argv[1]);

  printf("Benchmarking Traversal with N=%d particles...\n", N);

  // 1. Setup Data
  double *px = malloc(N * sizeof(double));
  double *py = malloc(N * sizeof(double));
  double *mass = malloc(N * sizeof(double));
  double *fx = malloc(N * sizeof(double));
  double *fy = malloc(N * sizeof(double));

  for (int i = 0; i < N; i++) {
    px[i] = (double)rand() / RAND_MAX;
    py[i] = (double)rand() / RAND_MAX;
    mass[i] = 1.0;
  }

  // 2. Build Tree
  NodeArena arena;
  init_arena(&arena, 3 * N);

  // Create Root
  TNode *root = create_new_TNode(&arena, -1, 0, 1, 0, 1);
  root->pos_x = px[0];
  root->pos_y = py[0];
  root->mass = mass[0];
  root->PID = 0;

  // Insert
  for (int i = 1; i < N; i++) {
    insert(&arena, root, px[i], py[i], mass[i], i);
  }
  printf("Tree Built. Nodes used: %zu\n", arena.used);

  // 3. Bench Recursive
  double t0 = get_time();
#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    double fxi = 0, fyi = 0;
    compute_force_recursive(px[i], py[i], mass[i], i, root, &fxi, &fyi, N, 0.5);
    fx[i] = fxi; // dummy write
  }
  double t_rec = get_time() - t0;
  printf("Recursive Time: %.4f s\n", t_rec);

  // 4. Bench Stackless
  t0 = get_time();
#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    double fxi = 0, fyi = 0;
    compute_force_stackless(px[i], py[i], mass[i], i, root, &fxi, &fyi, N, 0.5);
    fx[i] = fxi;
  }
  double t_stack = get_time() - t0;
  printf("Stackless Time: %.4f s\n", t_stack);

  printf("Speedup: %.2fx\n", t_rec / t_stack);

  // Cleanup
  free_arena(&arena);
  free(px);
  free(py);
  free(mass);
  free(fx);
  free(fy);
  return 0;
}
