#ifndef KERNELS_H
#define KERNELS_H

#include "../core/types.h"

// Kernel function pointer definition
typedef void (*ForceComputeKernel)(ParticleSystem *sys, KernelConfig *config);

// Kernel declarations
void compute_force_v1_naive(ParticleSystem *sys, KernelConfig *config);
void compute_force_v2_barnes_hut(ParticleSystem *sys, KernelConfig *config);
void compute_force_v3_arena(ParticleSystem *sys, KernelConfig *config);
void compute_force_v4_morton(ParticleSystem *sys, KernelConfig *config);
void compute_force_v5_parallel(ParticleSystem *sys, KernelConfig *config);

#endif
