#ifndef KMEANS_H
#define KMEANS_H

#include "ds.h"
#include "types.h"
#include <stdbool.h>

#define MAX_ITERATIONS 50

typedef struct {
  double ctr_x;
  double ctr_y;
  int count;
} CNode;

// Cluster particles into k groups.
bool kmeans(ParticleSystem *sys, int *clustersP, int *clusters_size, int k,
            int n_threads);

// Check centroid convergence.
bool converged(CNode *clusters, double *old_clusters_ctr_x,
               double *old_clusters_ctr_y, int iterations, int k);

// Update cluster centroids.
void getCentroids(ParticleSystem *sys, CNode *clusters, int *labels, int k,
                  int N);

// Assign each particle to the nearest cluster.
void assignLabels(ParticleSystem *sys, CNode *clusters, int *labels, int k,
                  int N, int n_threads);

#endif
