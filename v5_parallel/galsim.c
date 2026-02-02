#include "barnes_hut.h"
#include "kmeans.h"
#include "morton.h"

#define EPSILON_O 1e-3
#define CHUNK_SIZE 8
#define INTEGRATOR_ORDER 2 // 1: Euler, 2: Verlet, 4: RK4

#ifndef _OPENMP
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}
#endif

void get_accel(int N, double *px, double *py, double *mass,
               int *clusters, double *region, int *clusters_size, int k,
               double *fx, double *fy, double *mass_inver,
               double *ax, double *ay, int n_threads,
               double THETA_MAX, NodeArena *arena) {
  barnes_hut(px, py, mass, N, clusters, region, clusters_size, k, fx, fy,
             n_threads, THETA_MAX, arena);
  for (int i = 0; i < N; i++) {
    ax[i] = fx[i] * mass_inver[i];
    ay[i] = fy[i] * mass_inver[i];
  }
}

int main(int argc, char *argv[]) {
#ifdef _OPENMP
  double time_tol = omp_get_wtime();
#else
  double time_tol = get_wall_seconds();
#endif

  if (argc < 8) {
    printf("Usage: N filename nsteps delta_t n_threads theta_max k\n");
    return 1;
  }

  int N = atoi(argv[1]);
  char *filename = argv[2];
  int nsteps = atoi(argv[3]);
  double delta_t = atof(argv[4]);
  int n_threads = atoi(argv[5]);
  double THETA_MAX = atof(argv[6]);
  int k = atoi(argv[7]);

  /* 2D Particle Data */
  double *pos_x = (double *)malloc(N * sizeof(double));
  double *pos_y = (double *)malloc(N * sizeof(double));
  double *mass = (double *)malloc(N * sizeof(double));
  double *vx = (double *)malloc(N * sizeof(double));
  double *vy = (double *)malloc(N * sizeof(double));
  double *brightness = (double *)malloc(N * sizeof(double));

  double *fx = (double *)malloc(N * sizeof(double));
  double *fy = (double *)malloc(N * sizeof(double));
  double *acc_x = (double *)malloc(N * sizeof(double));
  double *acc_y = (double *)malloc(N * sizeof(double));
  double *mass_inver = (double *)malloc(N * sizeof(double));

  if (!pos_x || !pos_y || !mass || !vx || !vy || !brightness ||
      !fx || !fy || !acc_x || !acc_y || !mass_inver) {
    fprintf(stderr, "Memory allocation failed.\n");
    return 1;
  }

  /* Read data (Assuming 2D input) */
  FILE *data_file = fopen(filename, "rb");
  if (data_file == NULL) {
    printf("Error opening file!\n");
    return 1;
  }
  for (int i = 0; i < N; i++) {
    fread(&pos_x[i], sizeof(double), 1, data_file);
    fread(&pos_y[i], sizeof(double), 1, data_file);
    fread(&mass[i], sizeof(double), 1, data_file);
    mass_inver[i] = 1.0 / mass[i];
    fread(&vx[i], sizeof(double), 1, data_file);
    fread(&vy[i], sizeof(double), 1, data_file);
    fread(&brightness[i], sizeof(double), 1, data_file);
  }
  fclose(data_file);

  int *clusters = (int *)malloc(k * N * sizeof(int));
  int *clusters_size = (int *)malloc(k * sizeof(int));
  double region[4] = {pos_x[0], pos_x[0], pos_y[0], pos_y[0]};
  for (int i = 1; i < N; i++) {
    if (pos_x[i] < region[0]) region[0] = pos_x[i];
    if (pos_x[i] > region[1]) region[1] = pos_x[i];
    if (pos_y[i] < region[2]) region[2] = pos_y[i];
    if (pos_y[i] > region[3]) region[3] = pos_y[i];
  }
  // Add a small buffer to avoid particles on the boundary, and ensure a minimum width
  double dx = region[1] - region[0];
  double dy = region[3] - region[2];
  if (dx < 1e-6) dx = 1e-6;
  if (dy < 1e-6) dy = 1e-6;
  region[0] -= dx * 0.05; region[1] += dx * 0.05;
  region[2] -= dy * 0.05; region[3] += dy * 0.05;

  z_order_sort(pos_x, pos_y, mass, vx, vy, brightness, N, region[0],
               region[1], region[2], region[3]);

  if (k == 1) {
    clusters_size[0] = N;
    for (int i = 0; i < N; i++)
      clusters[i] = i;
  } else {
    kmeans(pos_x, pos_y, N, clusters, clusters_size, k, n_threads);
  }

  NodeArena arena;
  init_arena(&arena, 200 * N);

  /* Pre-allocate buffers for RK4 */
#if INTEGRATOR_ORDER == 4
  double *k1vx = (double *)malloc(N * sizeof(double));
  double *k1vy = (double *)malloc(N * sizeof(double));
  double *k1px = (double *)malloc(N * sizeof(double));
  double *k1py = (double *)malloc(N * sizeof(double));

  double *k2vx = (double *)malloc(N * sizeof(double));
  double *k2vy = (double *)malloc(N * sizeof(double));
  double *k2px = (double *)malloc(N * sizeof(double));
  double *k2py = (double *)malloc(N * sizeof(double));

  double *k3vx = (double *)malloc(N * sizeof(double));
  double *k3vy = (double *)malloc(N * sizeof(double));
  double *k3px = (double *)malloc(N * sizeof(double));
  double *k3py = (double *)malloc(N * sizeof(double));

  double *k4vx = (double *)malloc(N * sizeof(double));
  double *k4vy = (double *)malloc(N * sizeof(double));
  double *k4px = (double *)malloc(N * sizeof(double));
  double *k4py = (double *)malloc(N * sizeof(double));

  double *tmp_px = (double *)malloc(N * sizeof(double));
  double *tmp_py = (double *)malloc(N * sizeof(double));

  if (!k1vx || !k1vy || !k1px || !k1py ||
      !k2vx || !k2vy || !k2px || !k2py ||
      !k3vx || !k3vy || !k3px || !k3py ||
      !k4vx || !k4vy || !k4px || !k4py ||
      !tmp_px || !tmp_py) {
    fprintf(stderr, "RK4 Buffer allocation failed.\n");
    return 1;
  }
#endif

  /* Simulation Loop */
  FILE *movie_file = fopen("movie.gal", "wb");
  for (int step = 0; step < nsteps; step++) {
    // Save frame for animation
    if (step % 1 == 0) {
            for (int i = 0; i < N; i++) {
                fwrite(&pos_x[i], sizeof(double), 1, movie_file);
                fwrite(&pos_y[i], sizeof(double), 1, movie_file);
                fwrite(&mass[i], sizeof(double), 1, movie_file);
            }
        }

#if INTEGRATOR_ORDER == 4 // RK4
    /* Stage 1: k1 = f(t, y) */
    get_accel(N, pos_x, pos_y, mass, clusters, region, clusters_size, k,
              fx, fy, mass_inver, k1vx, k1vy, n_threads, THETA_MAX, &arena);
    for (int i = 0; i < N; i++) {
      k1px[i] = vx[i];
      k1py[i] = vy[i];
    }

    /* Stage 2: k2 = f(t + dt/2, y + dt/2 * k1) */
    for (int i = 0; i < N; i++) {
      tmp_px[i] = pos_x[i] + 0.5 * delta_t * k1px[i];
      tmp_py[i] = pos_y[i] + 0.5 * delta_t * k1py[i];
    }
    get_accel(N, tmp_px, tmp_py, mass, clusters, region, clusters_size, k,
              fx, fy, mass_inver, k2vx, k2vy, n_threads, THETA_MAX, &arena);
    for (int i = 0; i < N; i++) {
      k2px[i] = vx[i] + 0.5 * delta_t * k1vx[i];
      k2py[i] = vy[i] + 0.5 * delta_t * k1vy[i];
    }

    /* Stage 3: k3 = f(t + dt/2, y + dt/2 * k2) */
    for (int i = 0; i < N; i++) {
      tmp_px[i] = pos_x[i] + 0.5 * delta_t * k2px[i];
      tmp_py[i] = pos_y[i] + 0.5 * delta_t * k2py[i];
    }
    get_accel(N, tmp_px, tmp_py, mass, clusters, region, clusters_size, k,
              fx, fy, mass_inver, k3vx, k3vy, n_threads, THETA_MAX, &arena);
    for (int i = 0; i < N; i++) {
      k3px[i] = vx[i] + 0.5 * delta_t * k2vx[i];
      k3py[i] = vy[i] + 0.5 * delta_t * k2vy[i];
    }

    /* Stage 4: k4 = f(t + dt, y + dt * k3) */
    for (int i = 0; i < N; i++) {
      tmp_px[i] = pos_x[i] + delta_t * k3px[i];
      tmp_py[i] = pos_y[i] + delta_t * k3py[i];
    }
    get_accel(N, tmp_px, tmp_py, mass, clusters, region, clusters_size, k,
              fx, fy, mass_inver, k4vx, k4vy, n_threads, THETA_MAX, &arena);
    for (int i = 0; i < N; i++) {
      k4px[i] = vx[i] + delta_t * k3vx[i];
      k4py[i] = vy[i] + delta_t * k3vy[i];
    }

    /* Final Update: y = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4) */
    for (int i = 0; i < N; i++) {
      pos_x[i] += (delta_t / 6.0) * (k1px[i] + 2.0 * k2px[i] + 2.0 * k3px[i] + k4px[i]);
      pos_y[i] += (delta_t / 6.0) * (k1py[i] + 2.0 * k2py[i] + 2.0 * k3py[i] + k4py[i]);
      vx[i] += (delta_t / 6.0) * (k1vx[i] + 2.0 * k2vx[i] + 2.0 * k3vx[i] + k4vx[i]);
      vy[i] += (delta_t / 6.0) * (k1vy[i] + 2.0 * k2vy[i] + 2.0 * k3vy[i] + k4vy[i]);
    }
#endif

#if INTEGRATOR_ORDER == 2 // Velocity Verlet
    // 1. Half-step velocity update
    get_accel(N, pos_x, pos_y, mass, clusters, region, clusters_size, k,
              fx, fy, mass_inver, acc_x, acc_y, n_threads, THETA_MAX,
              &arena);
    for (int i = 0; i < N; i++) {
      vx[i] += 0.5 * delta_t * acc_x[i];
      vy[i] += 0.5 * delta_t * acc_y[i];
      pos_x[i] += delta_t * vx[i];
      pos_y[i] += delta_t * vy[i];
    }
    // 2. Full-step velocity update
    get_accel(N, pos_x, pos_y, mass, clusters, region, clusters_size, k,
              fx, fy, mass_inver, acc_x, acc_y, n_threads, THETA_MAX,
              &arena);
    for (int i = 0; i < N; i++) {
      vx[i] += 0.5 * delta_t * acc_x[i];
      vy[i] += 0.5 * delta_t * acc_y[i];
    }
#endif

#if INTEGRATOR_ORDER == 1 // Euler
    get_accel(N, pos_x, pos_y, mass, clusters, region, clusters_size, k,
              fx, fy, mass_inver, acc_x, acc_y, n_threads, THETA_MAX,
              &arena);
    for (int i = 0; i < N; i++) {
      pos_x[i] += delta_t * vx[i];
      pos_y[i] += delta_t * vy[i];
      vx[i] += delta_t * acc_x[i];
      vy[i] += delta_t * acc_y[i];
    }
#endif

    if (k > 1 && (step % 10 == 0)) {
      kmeans(pos_x, pos_y, N, clusters, clusters_size, k, n_threads);
    }
  }
  fclose(movie_file);

  /* Write Result */
  FILE *rfile = fopen("result.gal", "wb");
  for (int i = 0; i < N; i++) {
    fwrite(&pos_x[i], sizeof(double), 1, rfile);
    fwrite(&pos_y[i], sizeof(double), 1, rfile);
    fwrite(&mass[i], sizeof(double), 1, rfile);
    fwrite(&vx[i], sizeof(double), 1, rfile);
    fwrite(&vy[i], sizeof(double), 1, rfile);
    fwrite(&brightness[i], sizeof(double), 1, rfile);
  }
  fclose(rfile);

  printf("Simulation took %7.8f seconds.\n",
         (omp_get_wtime() - time_tol));

  free_arena(&arena);
  free(pos_x); free(pos_y);
  free(mass); free(mass_inver);
  free(vx); free(vy);
  free(brightness);
  free(fx); free(fy);
  free(acc_x); free(acc_y);
  free(clusters); free(clusters_size);

#if INTEGRATOR_ORDER == 4
  free(k1vx); free(k1vy);
  free(k1px); free(k1py);
  free(k2vx); free(k2vy);
  free(k2px); free(k2py);
  free(k3vx); free(k3vy);
  free(k3px); free(k3py);
  free(k4vx); free(k4vy);
  free(k4px); free(k4py);
  free(tmp_px); free(tmp_py);
#endif

  return 0;
}
