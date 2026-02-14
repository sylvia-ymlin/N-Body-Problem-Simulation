#include "barnes_hut.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int main(int argc, char *argv[]) {
  if (argc < 8) {
    printf("Usage: N filename nsteps delta_t n_threads theta_max k\n");
    return 1;
  }

  int N = atoi(argv[1]);
  char *filename = argv[2];
  int nsteps = atoi(argv[3]);
  double dt = atof(argv[4]);
  // n_threads and k are ignored in this version
  double theta = atof(argv[6]);

  double *px = malloc(N * sizeof(double));
  double *py = malloc(N * sizeof(double));
  double *mass = malloc(N * sizeof(double));
  double *vx = malloc(N * sizeof(double));
  double *vy = malloc(N * sizeof(double));
  double *fx = malloc(N * sizeof(double));
  double *fy = malloc(N * sizeof(double));

  FILE *f = fopen(filename, "rb");
  for (int i = 0; i < N; i++) {
    fread(&px[i], 8, 1, f);
    fread(&py[i], 8, 1, f);
    fread(&mass[i], 8, 1, f);
    fread(&vx[i], 8, 1, f);
    fread(&vy[i], 8, 1, f);

    fseek(f, 8, SEEK_CUR); // Skip brightness
  }
  fclose(f);

  NodeArena arena;
  init_arena(&arena, 1000 * N);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  FILE *movie_file = fopen("movie.gal", "wb");
  
  // Pre-compute initial forces (moved outside loop)
  barnes_hut(px, py, mass, N, fx, fy, theta, &arena);

  for (int step = 0; step < nsteps; step++) {
    if (step % 1 == 0) {
      for (int i = 0; i < N; i++) {
        fwrite(&px[i], 8, 1, movie_file);
        fwrite(&py[i], 8, 1, movie_file);
        fwrite(&mass[i], 8, 1, movie_file);
      }
    }
    
    // 1. First half-kick and Drift (Velocity Verlet)
    for (int i = 0; i < N; i++) {
      vx[i] += 0.5 * dt * fx[i] / mass[i];
      vy[i] += 0.5 * dt * fy[i] / mass[i];
      px[i] += dt * vx[i];
      py[i] += dt * vy[i];
    }

    // 2. Recompute Force
    barnes_hut(px, py, mass, N, fx, fy, theta, &arena);

    // 3. Final Velocity
    for (int i = 0; i < N; i++) {
      vx[i] += 0.5 * dt * fx[i] / mass[i];
      vy[i] += 0.5 * dt * fy[i] / mass[i];
    }
  }
  fclose(movie_file);

  gettimeofday(&end, NULL);
  double elapsed =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
  printf("Barnes-Hut (Arena) took %f seconds\n", elapsed);

  FILE *rfile = fopen("result.gal", "wb");
  if (rfile) {
    for (int i = 0; i < N; i++) {
      fwrite(&px[i], sizeof(double), 1, rfile);
      fwrite(&py[i], sizeof(double), 1, rfile);
      fwrite(&mass[i], sizeof(double), 1, rfile);
      fwrite(&vx[i], sizeof(double), 1, rfile);
      fwrite(&vy[i], sizeof(double), 1, rfile);
    }
    fclose(rfile);
  }

  free_arena(&arena);
  free(px);
  free(py);
  free(mass);
  free(vx);
  free(vy);
  free(fx);
  free(fy);
  return 0;
}
