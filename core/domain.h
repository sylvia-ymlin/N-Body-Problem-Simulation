#ifndef DOMAIN_H
#define DOMAIN_H

// Shared domain state for Barnes-Hut variants (one instance per kernel TU).
typedef struct {
    int initialized;
    int N;
    double x_min, x_max, y_min, y_max;
} Domain;

// Square root-domain helpers used by Barnes–Hut variants.
//
// The domain is derived from current particle positions as a square with side
// length d = max(dx, dy), then padded by (padding_frac * d) on all sides.
// Updates are expand-only (never shrink) to keep a stable global domain.

void domain_square_init_from_particles(const double* x, const double* y, int N,
                                       double padding_frac, double* x_min,
                                       double* x_max, double* y_min,
                                       double* y_max);

void domain_square_expand_if_needed(const double* x, const double* y, int N,
                                    double padding_frac, double* x_min,
                                    double* x_max, double* y_min,
                                    double* y_max);

#endif
