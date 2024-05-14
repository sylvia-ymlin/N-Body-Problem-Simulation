#include "core/io.h"
#include "core/time_utils.h"
#include "core/types.h"
#include "simulation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    int version;
    int N;
    const char* filename;
    int nsteps;
    double dt;
    int n_threads;
    double theta;
    int k_clusters;
} SimulationArgs;

static void print_usage(const char* program_name);
static int parse_args(int argc, char* argv[], SimulationArgs* args);
static ForceComputeKernel select_kernel(int version);
static void integrate_positions(ParticleSystem* sys, double dt,
                                int use_parallel);
static void integrate_velocities(ParticleSystem* sys, double dt,
                                 int use_parallel);
static void run_simulation(ParticleSystem* sys, const SimulationArgs* args,
                           ForceComputeKernel kernel, KernelConfig* config);

int main(int argc, char* argv[]) {
    SimulationArgs args;
    if (!parse_args(argc, argv, &args)) {
        print_usage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    omp_set_num_threads(args.n_threads);
#endif

    KernelConfig config = {
        args.theta,
        args.n_threads,
        args.k_clusters,
        0.0,
    };

    ForceComputeKernel kernel = select_kernel(args.version);
    if (!kernel) {
        fprintf(stderr, "Version %d not found\n", args.version);
        return 1;
    }

    ParticleSystem sys = io_read_particles(args.filename, args.N);
    run_simulation(&sys, &args, kernel, &config);

    char out_name[64];
    sprintf(out_name, "data/outputs/result_v%d.gal", args.version);
    io_write_result(out_name, &sys);
    io_free_particles(&sys);
    return 0;
}

static void print_usage(const char* program_name) {
    printf("Usage: %s <version> N <input.gal> nsteps dt n_threads theta k\n",
           program_name);
    printf("Versions:\n");
    printf("  1: Naive pairwise reference\n");
    printf("  2: Barnes-Hut tree\n");
    printf("  3: Barnes-Hut + arena allocator\n");
    printf("  4: Barnes-Hut + Morton ordering\n");
    printf("  5: V4 + OpenMP force loop\n");
}

static int parse_args(int argc, char* argv[], SimulationArgs* args) {
    if (argc < 9) {
        return 0;
    }

    args->version = atoi(argv[1]);
    args->N = atoi(argv[2]);
    args->filename = argv[3];
    args->nsteps = atoi(argv[4]);
    args->dt = atof(argv[5]);
    args->n_threads = atoi(argv[6]);
    args->theta = atof(argv[7]);
    args->k_clusters = atoi(argv[8]);
    return 1;
}

static ForceComputeKernel select_kernel(int version) {
    switch (version) {
    case 1:
        return compute_force_v1_naive;
    case 2:
        return compute_force_v2_barnes_hut;
    case 3:
        return compute_force_v3_arena;
    case 4:
        return compute_force_v4_morton;
    case 5:
        return compute_force_v5_parallel;
    default:
        return NULL;
    }
}

static void integrate_positions(ParticleSystem* sys, double dt,
                                int use_parallel) {
    int N = sys->N;
#pragma omp parallel for if (use_parallel)
    for (int i = 0; i < N; i++) {
        double m_inv = 1.0 / sys->mass[i];
        sys->vx[i] += 0.5 * dt * sys->fx[i] * m_inv;
        sys->vy[i] += 0.5 * dt * sys->fy[i] * m_inv;
        sys->pos_x[i] += dt * sys->vx[i];
        sys->pos_y[i] += dt * sys->vy[i];
    }
}

static void integrate_velocities(ParticleSystem* sys, double dt,
                                 int use_parallel) {
    int N = sys->N;
#pragma omp parallel for if (use_parallel)
    for (int i = 0; i < N; i++) {
        double m_inv = 1.0 / sys->mass[i];
        sys->vx[i] += 0.5 * dt * sys->fx[i] * m_inv;
        sys->vy[i] += 0.5 * dt * sys->fy[i] * m_inv;
    }
}

static void run_simulation(ParticleSystem* sys, const SimulationArgs* args,
                           ForceComputeKernel kernel, KernelConfig* config) {
    const int use_parallel = (args->version == 5);

    printf("Running v%d | N=%d | Steps=%d | Threads=%d\n", args->version,
           args->N, args->nsteps, args->n_threads);

    double start = sim_time_now();
    kernel(sys, config);

    for (int step = 0; step < args->nsteps; step++) {
        // Velocity Verlet: half kick, full drift, recompute forces, half kick.
        integrate_positions(sys, args->dt, use_parallel);

        config->current_time = (step + 1) * args->dt;
        kernel(sys, config);

        integrate_velocities(sys, args->dt, use_parallel);

        if (step % 50 == 0) {
            printf("Step %d/%d\r", step, args->nsteps);
        }
        fflush(stdout);
    }

    double end = sim_time_now();
    printf("\nSimulation Complete: %.2fs\n", end - start);
}
