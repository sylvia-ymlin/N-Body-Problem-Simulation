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

/* Function: kmeans
 * --------------------------------------------------
 * Clusters particles into k clusters.
 * Populates clustersP with particle indices grouped by cluster.
 * Populates clusters_size with size of each cluster.
 */
bool kmeans(ParticleSystem *sys, int *clustersP, int *clusters_size, int k,
            int n_threads);

/* Function: converged
 * --------------------------------------------------
 * Checks if centroids have stabilized.
 */
bool converged(CNode *clusters, double *old_clusters_ctr_x,
               double *old_clusters_ctr_y, int iterations, int k);

/* Function: getCentroids
 * --------------------------------------------------
 * Update cluster centroids based on current assignment.
 */
void getCentroids(ParticleSystem *sys, CNode *clusters, int *labels, int k,
                  int N);

/* Function: assignLabels
 * --------------------------------------------------
 * Assigns each particle to the nearest cluster.
 */
void assignLabels(ParticleSystem *sys, CNode *clusters, int *labels, int k,
                  int N, int n_threads);

#endif
