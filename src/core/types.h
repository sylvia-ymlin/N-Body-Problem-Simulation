#ifndef TYPES_H
#define TYPES_H

typedef struct {
    int N;
    double *pos_x;
    double *pos_y;
    double *mass;
    double *vx;
    double *vy;
    double *fx;     // Force x component
    double *fy;     // Force y component
} ParticleSystem;

// Kernel configuration context
typedef struct {
    double theta_max;
    int n_threads;
    int use_arena;
    int k_clusters;
} KernelConfig;

#endif
