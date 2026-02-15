#include "io.h"
#include <stdio.h>
#include <stdlib.h>

ParticleSystem io_read_particles(const char *filename, int N) {
  ParticleSystem sys;
  sys.N = N;
  sys.pos_x = malloc(N * sizeof(double));
  sys.pos_y = malloc(N * sizeof(double));
  sys.mass = malloc(N * sizeof(double));
  sys.vx = malloc(N * sizeof(double));
  sys.vy = malloc(N * sizeof(double));
  sys.fx = malloc(N * sizeof(double));
  sys.fy = malloc(N * sizeof(double));

  if (!sys.pos_x || !sys.pos_y || !sys.mass || !sys.vx || !sys.vy || !sys.fx ||
      !sys.fy) {
    fprintf(stderr, "Error: Memory allocation failed for %d particles.\n", N);
    exit(1);
  }

  FILE *f = fopen(filename, "rb");
  if (!f) {
    perror("Error opening input file");
    exit(1);
  }

  // Input format: [x, y, mass, vx, vy, brightness] (6 doubles) per particle
  for (int i = 0; i < N; i++) {
    double x, y, m, vx, vy;
    if (fread(&x, 8, 1, f) != 1 || fread(&y, 8, 1, f) != 1 ||
        fread(&m, 8, 1, f) != 1 || fread(&vx, 8, 1, f) != 1 ||
        fread(&vy, 8, 1, f) != 1) {
      fprintf(stderr, "Error: Unexpected end of file at particle %d\n", i);
      exit(1);
    }
    fseek(f, 8, SEEK_CUR); // Skip brightness

    sys.pos_x[i] = x;
    sys.pos_y[i] = y;
    sys.mass[i] = m;
    sys.vx[i] = vx;
    sys.vy[i] = vy;
    sys.fx[i] = 0.0;
    sys.fy[i] = 0.0;
  }

  fclose(f);
  return sys;
}

void io_write_frame(const char *filename, ParticleSystem *sys) {
  FILE *f = fopen(filename, "ab"); // Append mode
  if (!f)
    return;

  for (int i = 0; i < sys->N; i++) {
    fwrite(&sys->pos_x[i], 8, 1, f);
    fwrite(&sys->pos_y[i], 8, 1, f);
    fwrite(&sys->mass[i], 8, 1, f);
  }
  fclose(f);
}

void io_write_result(const char *filename, ParticleSystem *sys) {
  FILE *f = fopen(filename, "wb");
  if (!f)
    return;

  for (int i = 0; i < sys->N; i++) {
    fwrite(&sys->pos_x[i], 8, 1, f);
    fwrite(&sys->pos_y[i], 8, 1, f);
    fwrite(&sys->mass[i], 8, 1, f);
    fwrite(&sys->vx[i], 8, 1, f);
    fwrite(&sys->vy[i], 8, 1, f);
  }
  fclose(f);
}

void io_free_particles(ParticleSystem *sys) {
  free(sys->pos_x);
  free(sys->pos_y);
  free(sys->mass);
  free(sys->vx);
  free(sys->vy);
  free(sys->fx);
  free(sys->fy);
}
