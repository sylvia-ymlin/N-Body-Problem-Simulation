#include "barnes_hut.h"
#include "morton.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int main(int argc, char *argv[]) {
  if (argc < 4) {
    printf("Usage: N filename nsteps [theta]\n");
    return 1;
  }

  int N = atoi(argv[1]);
  char *filename = argv[2];
  int nsteps = atoi(argv[3]);
  double theta = (argc > 4) ? atof(argv[4]) : 0.5;
  double dt = 0.001;

  double *px = malloc(N * sizeof(double));
  double *py = malloc(N * sizeof(double));
  double *mass = malloc(N * sizeof(double));
  double *vx = malloc(N * sizeof(double));
  double *vy = malloc(N * sizeof(double));
  double *brightness = malloc(N * sizeof(double));
  double *fx = malloc(N * sizeof(double));
  double *fy = malloc(N * sizeof(double));

  FILE *f = fopen(filename, "rb");
  if (!f) {
    printf("Error opening file!\n");
    return 1;
  }
  for (int i = 0; i < N; i++) {
    fread(&px[i], 8, 1, f); fread(&py[i], 8, 1, f);
    fread(&mass[i], 8, 1, f); fread(&vx[i], 8, 1, f);
    fread(&vy[i], 8, 1, f); fread(&brightness[i], 8, 1, f);
  }
  fclose(f);

  // Initial Morton sorting for spatial locality
  double x_min = px[0], x_max = px[0], y_min = py[0], y_max = py[0];
  for (int i = 1; i < N; i++) {
    if (px[i] < x_min) x_min = px[i];
    if (px[i] > x_max) x_max = px[i];
    if (py[i] < y_min) y_min = py[i];
    if (py[i] > y_max) y_max = py[i];
  }
  z_order_sort(px, py, mass, vx, vy, brightness, N, x_min, x_max, y_min, y_max);

  NodeArena arena;
  init_arena(&arena, 200 * N);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  for (int step = 0; step < nsteps; step++) {
    barnes_hut(px, py, mass, N, fx, fy, theta, &arena);
    for (int i = 0; i < N; i++) {
      vx[i] += dt * fx[i] / mass[i];
      vy[i] += dt * fy[i] / mass[i];
      px[i] += dt * vx[i];
      py[i] += dt * vy[i];
    }
    // Optionally re-sort every few steps to maintain locality
    if (step % 10 == 0 && step > 0) {
      z_order_sort(px, py, mass, vx, vy, brightness, N, x_min, x_max, y_min, y_max);
    }
  }

  gettimeofday(&end, NULL);
  double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
  printf("Barnes-Hut (Arena + Morton) took %f seconds\n", elapsed);

  free_arena(&arena);
  free(px); free(py); free(mass); free(vx); free(vy); free(brightness); free(fx); free(fy);
  return 0;
}
