#ifndef SIMULATION_H
#define SIMULATION_H

#include "src/types.h"

// Core simulation interface
typedef void (*ForceComputeKernel)(ParticleSystem *sys, KernelConfig *config);

// Solver implementations (v1-v5)
void compute_force_v1_naive(ParticleSystem *sys, KernelConfig *config);
void compute_force_v2_barnes_hut(ParticleSystem *sys, KernelConfig *config);
void compute_force_v3_arena(ParticleSystem *sys, KernelConfig *config);
void compute_force_v4_morton(ParticleSystem *sys, KernelConfig *config);
void compute_force_v5_parallel(ParticleSystem *sys, KernelConfig *config);

#define EPSILON 1e-3
#include <stddef.h>
#include <sys/time.h>
static inline double get_sim_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}
#endif
