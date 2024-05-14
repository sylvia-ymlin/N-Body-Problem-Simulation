#ifndef TYPES_H
#define TYPES_H

typedef struct {
    int N;
    double* pos_x;
    double* pos_y;
    double* mass;
    double* vx;
    double* vy;
    double* fx;
    double* fy;
} ParticleSystem;

typedef struct {
    double theta_max;
    int    n_threads;
    int    k_clusters;  /* 0 = use Morton ordering, >0 = use k-means clustering */
    double current_time;
} KernelConfig;

#endif
