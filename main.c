#include "io.h"
#include "time_utils.h"
#include "types.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void compute_force_naive(ParticleSystem* sys, KernelConfig* config);
void compute_force_barnes_hut(ParticleSystem* sys, KernelConfig* config);

static void integrate_positions(ParticleSystem* sys, double dt);
static void integrate_velocities(ParticleSystem* sys, double dt);

static const int    DEFAULT_STEPS   = 200;
static const double DEFAULT_DT      = 1e-5;
static const int    DEFAULT_THREADS = 8;
static const double DEFAULT_THETA   = 0.5;
static const int    DEFAULT_K       = 0;

int main(int argc, char* argv[]) {
    if (argc != 4 && argc != 5 && argc != 6 && argc != 9) {
        fprintf(stderr, "Usage: %s <version> N <input.gal> [nsteps] [n_threads]\n", argv[0]);
        fprintf(stderr, "   or: %s <version> N <input.gal> nsteps dt n_threads theta k\n", argv[0]);
        fprintf(stderr, "Versions: 1=Naive  2=Barnes-Hut\n");
        fprintf(stderr, "k: locality strategy for version 2 (0=Morton, >0=k-means with k clusters)\n");
        return 1;
    }

    char* end = NULL;
    int version_id = (int)strtol(argv[1], &end, 10);
    if (end == argv[1] || *end != '\0' || version_id < 1 || version_id > 2) {
        fprintf(stderr, "Version must be 1 or 2.\n");
        return 1;
    }

    int N = (int)strtol(argv[2], &end, 10);
    if (end == argv[2] || *end != '\0' || N <= 0) return 1;

    const char* filename = argv[3];
    int    nsteps     = DEFAULT_STEPS;
    double dt         = DEFAULT_DT;
    int    n_threads  = DEFAULT_THREADS;
    double theta      = DEFAULT_THETA;
    int    k_clusters = DEFAULT_K;

    if (argc >= 5) {
        nsteps = (int)strtol(argv[4], &end, 10);
        if (end == argv[4] || *end != '\0' || nsteps <= 0) return 1;
    }
    if (argc >= 6) {
        n_threads = (int)strtol(argv[5], &end, 10);
        if (end == argv[5] || *end != '\0' || n_threads <= 0) return 1;
    }
    if (argc == 9) {
        dt = strtod(argv[5], &end);
        if (end == argv[5] || *end != '\0' || dt <= 0.0) return 1;
        n_threads = (int)strtol(argv[6], &end, 10);
        if (end == argv[6] || *end != '\0' || n_threads <= 0) return 1;
        theta = strtod(argv[7], &end);
        if (end == argv[7] || *end != '\0' || theta < 0.0) return 1;
        k_clusters = (int)strtol(argv[8], &end, 10);
        if (end == argv[8] || *end != '\0' || k_clusters < 0) return 1;
    }

    if (version_id == 1 && k_clusters != 0) {
        fprintf(stderr, "k-means is only supported for version 2.\n");
        return 1;
    }

#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    KernelConfig config = { theta, n_threads, k_clusters, 0.0 };
    ParticleSystem sys = io_read_particles(filename, N);

    printf("Running v%d | N=%d | Steps=%d | Threads=%d\n",
           version_id, sys.N, nsteps, n_threads);
    printf("dt=%.1e | theta=%.2f | k=%d\n", dt, theta, k_clusters);

    /* Initial force computation */
    if (version_id == 1) compute_force_naive(&sys, &config);
    else                 compute_force_barnes_hut(&sys, &config);

    double t_start = sim_time_now();

    for (int step = 0; step < nsteps; step++) {
        /* Velocity Verlet: half kick, full drift, recompute forces, half kick */
        integrate_positions(&sys, dt);

        config.current_time = (step + 1) * dt;

        if (version_id == 1) compute_force_naive(&sys, &config);
        else                 compute_force_barnes_hut(&sys, &config);

        integrate_velocities(&sys, dt);

        if (step % 50 == 0)
            printf("Step %d/%d\r", step, nsteps);
        fflush(stdout);
    }

    printf("\nSimulation Complete: %.2fs\n", sim_time_now() - t_start);

    char out_name[64];
    const char* label = (version_id == 1) ? "naive" : "barnes_hut";
    snprintf(out_name, sizeof(out_name), "data/outputs/result_%s.gal", label);
    io_write_result(out_name, &sys);
    io_free_particles(&sys);
    return 0;
}

/** Velocity Verlet half-kick + full drift
 * Updates velocities by half a step then advances positions.
 * ----------------------------------------------------------------- */
static void integrate_positions(ParticleSystem* sys, double dt) {
    int N = sys->N;
    for (int i = 0; i < N; i++) {
        double m_inv = 1.0 / sys->mass[i];
        sys->vx[i] += 0.5 * dt * sys->fx[i] * m_inv;
        sys->vy[i] += 0.5 * dt * sys->fy[i] * m_inv;
        sys->pos_x[i] += dt * sys->vx[i];
        sys->pos_y[i] += dt * sys->vy[i];
    }
}

/** Velocity Verlet second half-kick
 * Completes the velocity update after new forces are computed.
 * ----------------------------------------------------------------- */
static void integrate_velocities(ParticleSystem* sys, double dt) {
    int N = sys->N;
    for (int i = 0; i < N; i++) {
        double m_inv = 1.0 / sys->mass[i];
        sys->vx[i] += 0.5 * dt * sys->fx[i] * m_inv;
        sys->vy[i] += 0.5 * dt * sys->fy[i] * m_inv;
    }
}
