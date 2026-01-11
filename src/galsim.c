#include "barnes_hut.h"
#include "kmeans.h"
#include "morton.h"

#define EPSILON_O 1e-3
#define CHUNK_SIZE 8
#define TWO_ORDER 1

#ifndef _OPENMP
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}
#endif

int main(int argc, char *argv[]) {
#ifdef _OPENMP
  double time_tol = omp_get_wtime();
#else
  double time_tol = get_wall_seconds();
#endif

  if (argc != 8) {
    printf("You should enter the following parameters in order:\n");
    printf("N filname nsteps delta_t n_threads theta_max k\n");
    return 1;
  }

  /* Accept the parameters */
  int N = atoi(argv[1]);
  char *filename = argv[2];
  int nsteps = atoi(argv[3]);
  double delta_t = atof(argv[4]);
  int n_threads = atoi(argv[5]);
  double THETA_MAX = atof(argv[6]);
  int k = atoi(argv[7]);

  /* Variables for the simulation
   * -------------------------------------------------------------------------------------*/
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

  // Check for allocation failures
  if (!pos_x || !pos_y || !mass || !vx || !vy || !brightness || !fx || !fy ||
      !acc_x || !acc_y || !mass_inver) {
    fprintf(stderr, "Memory allocation failed for particle arrays.\n");
    return 1;
  }

  double new_acc_x, new_acc_y;

  /* Read data from the file
   * -------------------------------------------------------------------------------------*/
  FILE *data_file = fopen(filename, "rb");
  if (data_file == NULL) {
    printf("Error opening file!\n");
    // Free memory before exit
    free(pos_x);
    free(pos_y);
    free(mass);
    free(vx);
    free(vy);
    free(brightness);
    free(fx);
    free(fy);
    free(acc_x);
    free(acc_y);
    free(mass_inver);
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

  /* Initialize the clusters
   * -------------------------------------------------------------------------------------*/
  // Allocate flat array on heap: size k * N
  int *clusters = (int *)malloc(k * N * sizeof(int));
  int *clusters_size = (int *)malloc(k * sizeof(int));

  if (!clusters || !clusters_size) {
    fprintf(stderr, "Memory allocation failed for clusters.\n");
    return 1;
  }

  double region[4] = {-100000, 100000, -100000,
                      100000}; // left, right, down, up

  /* Spatial Sorting (Z-Order) for Cache Locality
     This sorts the particle arrays in-place based on Morton codes to improve
     memory access patterns during tree traversal. */
  z_order_sort(pos_x, pos_y, mass, vx, vy, brightness, N, region[0], region[1],
               region[2], region[3]);

  if (k == 1) { /* No clustering */
    clusters_size[0] = N;
    for (int i = 0; i < N; i++) {
      // Flat array access: row 0, col i -> 0*N + i = i
      clusters[i] = i;
    }
  } else {
    kmeans(pos_x, pos_y, N, clusters, clusters_size, k, n_threads);
  }
  /* Initialize Memory Arena for Barnes-Hut Tree */
  NodeArena arena;
  // Capacity guess: 2*N is usually enough for Barnes-Hut.
  // If not, the arena allocator will exit(1).
  init_arena(&arena, 3 * N);

/* Velocity Verlet: initial acceleration
 * -------------------------------------------------------------------------------------*/
#if TWO_ORDER
  barnes_hut(pos_x, pos_y, mass, N, clusters, region, clusters_size, k, fx, fy,
             n_threads, THETA_MAX, &arena);
  for (int i = 0; i < N; i++) {
    acc_x[i] = fx[i] * mass_inver[i];
    acc_y[i] = fy[i] * mass_inver[i];
  }
#endif

  /* Simulation
   * -------------------------------------------------------------------------------------*/
  for (int step = 0; step < nsteps; step++) {
#if TWO_ORDER
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
#endif
    // - update the position first
    for (int i = 0; i < N; i++) {
      pos_x[i] += delta_t * vx[i] + 0.5 * delta_t * delta_t * acc_x[i];
      pos_y[i] += delta_t * vy[i] + 0.5 * delta_t * delta_t * acc_y[i];
      if (pos_x[i] < region[0] || pos_x[i] > region[1] ||
          pos_y[i] < region[2] || pos_y[i] > region[3]) {
        printf("At least one particle is out of the region, and the simulation "
               "has been terminated.\n");
        free_arena(&arena);
        exit(0);
      }
    }
#endif

    /* update clusters each 0.001 seconds */
    k > 1 && step != 0 && ((int)(step * delta_t * 100000) % 10 == 0) &&
        kmeans(pos_x, pos_y, N, clusters, clusters_size, k, n_threads);

    /* forces calculation */
    barnes_hut(pos_x, pos_y, mass, N, clusters, region, clusters_size, k, fx,
               fy, n_threads, THETA_MAX, &arena);

/* update the velocity and acceleration */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
#endif
    for (int i = 0; i < N; i++) {
#if TWO_ORDER
      new_acc_x = fx[i] * mass_inver[i];
      new_acc_y = fy[i] * mass_inver[i];
      vx[i] += 0.5 * delta_t * (new_acc_x + acc_x[i]);
      vy[i] += 0.5 * delta_t * (new_acc_y + acc_y[i]);
      acc_x[i] = new_acc_x;
      acc_y[i] = new_acc_y;
#else
      vx[i] += delta_t * fx[i] * mass_inver[i];
      vy[i] += delta_t * fy[i] * mass_inver[i];
      pos_x[i] += delta_t * vx[i];
      pos_y[i] += delta_t * vy[i];
      if (pos_x[i] < region[0] || pos_x[i] > region[1] ||
          pos_y[i] < region[2] || pos_y[i] > region[3]) {
        printf("At least one particle is out of the region, and the simulation "
               "has been terminated.\n");
        free_arena(&arena);
        exit(0);
      }
#endif
    }
  }

  /* write data into the file
   * -------------------------------------------------------------------------------------*/
  FILE *rfile = fopen("result.gal", "w");
  if (rfile == NULL) {
    printf("Error opening file!\n");
    free_arena(&arena);
    return 1;
  }
  for (int i = 0; i < N; i++) {
    fwrite(&pos_x[i], sizeof(double), 1, rfile);
    fwrite(&pos_y[i], sizeof(double), 1, rfile);
    fwrite(&mass[i], sizeof(double), 1, rfile);
    fwrite(&vx[i], sizeof(double), 1, rfile);
    fwrite(&vy[i], sizeof(double), 1, rfile);
    fwrite(&brightness[i], sizeof(double), 1, rfile);
  }
  fclose(rfile);

#ifdef _OPENMP
  printf("f_std tests took %7.8f wall seconds.\n", omp_get_wtime() - time_tol);
#else
  printf("f_std tests took %7.8f wall seconds.\n",
         get_wall_seconds() - time_tol);
#endif

  // Clean up memory
  free_arena(&arena);
  free(pos_x);

  free(pos_y);
  free(mass);
  free(vx);
  free(vy);
  free(brightness);
  free(fx);
  free(fy);
  free(acc_x);
  free(acc_y);
  free(mass_inver);
  free(clusters);
  free(clusters_size);

  return 0;
}
