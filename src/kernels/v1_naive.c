#include "kernels.h"
#include <math.h>
#include <stdlib.h>

void compute_force_v1_naive(ParticleSystem *sys, KernelConfig *config) {
  int N = sys->N;
  double eps = 1e-3;
  double G = 100.0 / N;

  // O(N^2) naive calculation
  for (int i = 0; i < N; i++) {
    double fx = 0.0;
    double fy = 0.0;
    for (int j = 0; j < N; j++) {
      if (i == j)
        continue;
      double dx = sys->pos_x[j] - sys->pos_x[i];
      double dy = sys->pos_y[j] - sys->pos_y[i];
      double distSq = dx * dx + dy * dy + eps * eps;
      double invDist = 1.0 / sqrt(distSq);
      double invDist3 = invDist * invDist * invDist;
      fx += G * sys->mass[i] * sys->mass[j] * dx * invDist3;
      fy += G * sys->mass[i] * sys->mass[j] * dy * invDist3;
    }
    sys->fx[i] = fx;
    sys->fy[i] = fy;
  }
}
