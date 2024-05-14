#include "simulation.h"
#include <math.h>

/* v1: Brute-force baseline, O(N^2) complexity. */
void compute_force_v1_naive(ParticleSystem* sys, KernelConfig* config) {
    (void)config;
    const int N = sys->N;
    const double G = G_FACTOR / N;

    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* m = sys->mass;
    double* fx_out = sys->fx;
    double* fy_out = sys->fy;

    for (int i = 0; i < N; i++) {
        double fx = 0.0;
        double fy = 0.0;

        for (int j = 0; j < N; j++) {
            if (i == j)
                continue;

            const double dx = x[j] - x[i];
            const double dy = y[j] - y[i];
            const double r = sqrt(dx * dx + dy * dy);
            const double denom = r + EPSILON;

            const double f = G * m[i] * m[j] / (denom * denom * denom);
            fx += f * dx;
            fy += f * dy;
        }
        fx_out[i] = fx;
        fy_out[i] = fy;
    }
}
