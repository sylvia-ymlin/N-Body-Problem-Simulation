#include "barnes_hut.h"
#include "kmeans.h"
#include "morton.h"

// Helper: Compute forces and update accelerations
// Inlined for clarity in the main loop
void compute_accelerations(int N, double *pos_x, double *pos_y, double *mass,
                           double *mass_inver, int *clusters,
                           int *clusters_size, int k, int n_threads,
                           double theta_max, NodeArena *arena, double *fx,
                           double *fy, double *ax, double *ay, double *region) {
  barnes_hut(pos_x, pos_y, mass, N, clusters, region, clusters_size, k, fx, fy,
             n_threads, theta_max, arena);
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
  // brightness removed (unused visualization artifact)

  double *fx = (double *)malloc(N * sizeof(double));
  double *fy = (double *)malloc(N * sizeof(double));
  double *acc_x = (double *)malloc(N * sizeof(double));
  double *acc_y = (double *)malloc(N * sizeof(double));
  double *mass_inver = (double *)malloc(N * sizeof(double));

  if (!pos_x || !pos_y || !mass || !vx || !vy || !fx || !fy || !acc_x ||
      !acc_y || !mass_inver) {
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
    fread(&vx[i], sizeof(double), 1, data_file);
    fread(&vy[i], sizeof(double), 1, data_file);
    fseek(data_file, sizeof(double), SEEK_CUR); // Skip brightness
  }
  fclose(data_file);

  int *clusters = (int *)malloc(N * sizeof(int));
  int *clusters_size = (int *)malloc(k * sizeof(int));
  double region[4] = {pos_x[0], pos_x[0], pos_y[0], pos_y[0]};
  for (int i = 1; i < N; i++) {
    if (pos_x[i] < region[0])
      region[0] = pos_x[i];
    if (pos_x[i] > region[1])
      region[1] = pos_x[i];
    if (pos_y[i] < region[2])
      region[2] = pos_y[i];
    if (pos_y[i] > region[3])
      region[3] = pos_y[i];
  }
  // Add a small buffer to avoid particles on the boundary, and ensure a minimum
  // width
  double dx = region[1] - region[0];
  double dy = region[3] - region[2];
  if (dx < 1e-6)
    dx = 1e-6;
  if (dy < 1e-6)
    dy = 1e-6;
  region[0] -= dx * 0.05;
  region[1] += dx * 0.05;
  region[2] -= dy * 0.05;
  region[3] += dy * 0.05;

  z_order_sort(pos_x, pos_y, mass, vx, vy, N, region[0], region[1], region[2],
               region[3]);

  for (int i = 0; i < N; i++) {
    mass_inver[i] = 1.0 / mass[i];
  }
  // WARN: mass_inver is not re-sorted if particles are re-ordered later!
  // It assumes 1:1 mapping with mass[] array.

  if (k == 1) {
    clusters_size[0] = N;
    for (int i = 0; i < N; i++)
      clusters[i] = i;
  } else {
    kmeans(pos_x, pos_y, N, clusters, clusters_size, k, n_threads);
  }

  NodeArena arena;
  init_arena(&arena, 1000 * N);

  /* Pre-allocate buffers for RK4 */
  // RK4 buffers removed. We use Velocity Verlet for symplectic stability.

  /* Simulation Loop */
  FILE *movie_file = fopen("movie.gal", "wb");
  for (int step = 0; step < nsteps; step++) {
    // Save frame for animation (sampled)
    if (step % 10 == 0) {
      for (int i = 0; i < N; i++) {
        fwrite(&pos_x[i], sizeof(double), 1, movie_file);
        fwrite(&pos_y[i], sizeof(double), 1, movie_file);
        fwrite(&mass[i], sizeof(double), 1, movie_file);
      }
    }

    /* Velocity Verlet Integration */
    // 1. First half-kick: v += 0.5 * a * dt
    // 2. Drift: x += v * dt
    compute_accelerations(N, pos_x, pos_y, mass, mass_inver, clusters,
                          clusters_size, k, n_threads, THETA_MAX, &arena, fx,
                          fy, acc_x, acc_y, region);

    for (int i = 0; i < N; i++) {
      vx[i] += 0.5 * delta_t * acc_x[i];
      vy[i] += 0.5 * delta_t * acc_y[i];
      pos_x[i] += delta_t * vx[i];
      pos_y[i] += delta_t * vy[i];
    }

    // 3. Recompute forces at new positions
    compute_accelerations(N, pos_x, pos_y, mass, mass_inver, clusters,
                          clusters_size, k, n_threads, THETA_MAX, &arena, fx,
                          fy, acc_x, acc_y, region);

    // 4. Second half-kick: v += 0.5 * a * dt
    for (int i = 0; i < N; i++) {
      vx[i] += 0.5 * delta_t * acc_x[i];
      vy[i] += 0.5 * delta_t * acc_y[i];
    }

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
    // fwrite(&brightness[i], sizeof(double), 1, rfile); // Removed
  }
  fclose(rfile);

#ifdef _OPENMP
  printf("Simulation took %7.8f seconds.\n", (omp_get_wtime() - time_tol));
#else
  printf("Simulation took %7.8f seconds.\n", (get_wall_seconds() - time_tol));
#endif

  free_arena(&arena);
  free(pos_x);
  free(pos_y);
  free(mass);
  free(mass_inver);
  free(vx);
  free(vy);
  // free(brightness); // Removed
  free(fx);
  free(fy);
  free(acc_x);
  free(acc_y);
  free(clusters);
  free(clusters_size);

  return 0;
}
