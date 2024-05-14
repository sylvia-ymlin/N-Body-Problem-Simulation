#include "domain.h"

void domain_square_init_from_particles(const double* x, const double* y, int N,
                                       double padding_frac, double* x_min,
                                       double* x_max, double* y_min,
                                       double* y_max) {
    double nx_min = x[0], nx_max = x[0];
    double ny_min = y[0], ny_max = y[0];
    for (int i = 1; i < N; i++) {
        if (x[i] < nx_min)
            nx_min = x[i];
        if (x[i] > nx_max)
            nx_max = x[i];
        if (y[i] < ny_min)
            ny_min = y[i];
        if (y[i] > ny_max)
            ny_max = y[i];
    }

    double dx = nx_max - nx_min;
    double dy = ny_max - ny_min;
    double d = (dx > dy) ? dx : dy;
    if (d <= 0.0)
        d = 1e-6;

    double sx_min = nx_min;
    double sx_max = nx_min + d;
    double sy_min = ny_min;
    double sy_max = ny_min + d;

    sx_min -= d * padding_frac;
    sx_max += d * padding_frac;
    sy_min -= d * padding_frac;
    sy_max += d * padding_frac;

    *x_min = sx_min;
    *x_max = sx_max;
    *y_min = sy_min;
    *y_max = sy_max;
}

void domain_square_expand_if_needed(const double* x, const double* y, int N,
                                    double padding_frac, double* x_min,
                                    double* x_max, double* y_min,
                                    double* y_max) {
    int needs_update = 0;
    for (int i = 0; i < N; i++) {
        if (x[i] < *x_min || x[i] > *x_max || y[i] < *y_min || y[i] > *y_max) {
            needs_update = 1;
            break;
        }
    }
    if (!needs_update)
        return;

    double sx_min, sx_max, sy_min, sy_max;
    domain_square_init_from_particles(x, y, N, padding_frac, &sx_min, &sx_max,
                                      &sy_min, &sy_max);

    if (sx_min < *x_min)
        *x_min = sx_min;
    if (sx_max > *x_max)
        *x_max = sx_max;
    if (sy_min < *y_min)
        *y_min = sy_min;
    if (sy_max > *y_max)
        *y_max = sy_max;
}
