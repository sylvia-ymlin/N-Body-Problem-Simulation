#include "barnes_hut.h"
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
  double *fx = malloc(N * sizeof(double));
  double *fy = malloc(N * sizeof(double));

  FILE *f = fopen(filename, "rb");
  for (int i = 0; i < N; i++) {
    fread(&px[i], 8, 1, f); fread(&py[i], 8, 1, f);
    fread(&mass[i], 8, 1, f); fread(&vx[i], 8, 1, f);
    fread(&vy[i], 8, 1, f); double b; fread(&b, 8, 1, f);
  }
  fclose(f);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  FILE *movie_file = fopen("movie.gal", "wb");
  for (int step = 0; step < nsteps; step++) {
    if (step % 1 == 0) {
      for (int i = 0; i < N; i++) {
        fwrite(&px[i], 8, 1, movie_file);
        fwrite(&py[i], 8, 1, movie_file);
        fwrite(&mass[i], 8, 1, movie_file);
      }
    }
    barnes_hut(px, py, mass, N, fx, fy, theta);
    for (int i = 0; i < N; i++) {
      vx[i] += dt * fx[i] / mass[i];
      vy[i] += dt * fy[i] / mass[i];
      px[i] += dt * vx[i];
      py[i] += dt * vy[i];
    }
  }
  fclose(movie_file);

  gettimeofday(&end, NULL);
  double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
  printf("Barnes-Hut (malloc) took %f seconds\n", elapsed);

  free(px); free(py); free(mass); free(vx); free(vy); free(fx); free(fy);
  return 0;
}
