#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
  double x, y, mass, vx, vy;
} Particle;

int main(int argc, char *argv[]) {
  if (argc < 8) {
    printf("Usage: N filename nsteps delta_t n_threads theta_max k\n");
    return 1;
  }

  int N = atoi(argv[1]);
  char *filename = argv[2];
  int nsteps = atoi(argv[3]);
  double dt = atof(argv[4]);
  // Note: n_threads, theta_max, and k are ignored in this naive version
  // but kept for interface consistency.

  double G = 100.0 / N;
  double eps = 1e-3;

  Particle *p = malloc(N * sizeof(Particle));
  FILE *f = fopen(filename, "rb");
  if (f == NULL) {
    printf("Error opening file!\n");
    return 1;
  }
  for (int i = 0; i < N; i++) {
    fread(&p[i].x, 8, 1, f);
    fread(&p[i].y, 8, 1, f);
    fread(&p[i].mass, 8, 1, f);
    fread(&p[i].vx, 8, 1, f);
    fread(&p[i].vy, 8, 1, f);
    fseek(f, 8, SEEK_CUR); // Skip brightness
  }
  fclose(f);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  FILE *movie_file = fopen("movie.gal", "wb");
  /* Simulation Loop */
  for (int s = 0; s < nsteps; s++) {
    // Save frame for animation
    if (s % 10 == 0) {
      for (int i = 0; i < N; i++) {
        fwrite(&p[i].x, 8, 1, movie_file);
        fwrite(&p[i].y, 8, 1, movie_file);
        fwrite(&p[i].mass, 8, 1, movie_file);
      }
    }

    // Velocity Verlet - Step 1: Half-step velocity and update positions
    // Since we need accelerations twice per step, let's keep it simple
    double *ax = malloc(N * sizeof(double));
    double *ay = malloc(N * sizeof(double));

    // 1a. Compute initial accelerations
    for (int i = 0; i < N; i++) {
      double fx = 0, fy = 0;
      for (int j = 0; j < N; j++) {
        if (i == j)
          continue;
        double dx = p[j].x - p[i].x;
        double dy = p[j].y - p[i].y;
        double distSq = dx * dx + dy * dy + eps * eps;
        double invDist = 1.0 / sqrt(distSq);
        double invDist3 = invDist * invDist * invDist;
        fx += G * p[i].mass * p[j].mass * dx * invDist3;
        fy += G * p[i].mass * p[j].mass * dy * invDist3;
      }
      ax[i] = fx / p[i].mass;
      ay[i] = fy / p[i].mass;
    }

    // 1b. Update particles
    for (int i = 0; i < N; i++) {
      p[i].vx += 0.5 * dt * ax[i];
      p[i].vy += 0.5 * dt * ay[i];
      p[i].x += dt * p[i].vx;
      p[i].y += dt * p[i].vy;
    }

    // 2a. Compute new accelerations at updated positions
    for (int i = 0; i < N; i++) {
      double fx = 0, fy = 0;
      for (int j = 0; j < N; j++) {
        if (i == j)
          continue;
        double dx = p[j].x - p[i].x;
        double dy = p[j].y - p[i].y;
        double distSq = dx * dx + dy * dy + eps * eps;
        double invDist = 1.0 / sqrt(distSq);
        double invDist3 = invDist * invDist * invDist;
        fx += G * p[i].mass * p[j].mass * dx * invDist3;
        fy += G * p[i].mass * p[j].mass * dy * invDist3;
      }
      ax[i] = fx / p[i].mass;
      ay[i] = fy / p[i].mass;
    }

    // 2b. Final half-step velocity update
    for (int i = 0; i < N; i++) {
      p[i].vx += 0.5 * dt * ax[i];
      p[i].vy += 0.5 * dt * ay[i];
    }

    free(ax);
    free(ay);
  }
  fclose(movie_file);

  gettimeofday(&end, NULL);
  double elapsed =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
  printf("Naive Simulation took %f seconds\n", elapsed);

  /* Write Result */
  FILE *rfile = fopen("result.gal", "wb");
  if (rfile) {
    for (int i = 0; i < N; i++) {
      fwrite(&p[i].x, 8, 1, rfile);
      fwrite(&p[i].y, 8, 1, rfile);
      fwrite(&p[i].mass, 8, 1, rfile);
      fwrite(&p[i].vx, 8, 1, rfile);
      fwrite(&p[i].vy, 8, 1, rfile);
    }
    fclose(rfile);
  }

  free(p);
  return 0;
}
