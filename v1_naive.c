#include "simulation.h"
#include <math.h>
#include <stdlib.h>

// v1: Brute-force baseline, O(N^2).
void compute_force_v1_naive(ParticleSystem *sys, KernelConfig *config) {
  int N = sys->N;
  double G = 100.0 / N; // Gravitational constant scaled by system size

  for (int i = 0; i < N; i++) {
    double fx = 0.0;
    double fy = 0.0;
    for (int j = 0; j < N; j++) {
      if (i == j)
        continue;

      double dx = sys->pos_x[j] - sys->pos_x[i];
      double dy = sys->pos_y[j] - sys->pos_y[i];
      double distSq = dx * dx + dy * dy + EPSILON * EPSILON;
      double invDist = 1.0 / sqrt(distSq);
      double invDist3 = invDist * invDist * invDist;

      double f = G * sys->mass[i] * sys->mass[j] * invDist3;
      fx += f * dx;
      fy += f * dy;
    }
    sys->fx[i] = fx;
    sys->fy[i] = fy;
  }
}
