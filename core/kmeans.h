#ifndef KMEANS_H
#define KMEANS_H

#include "types.h"
#include <stdbool.h>

// Cluster particles into k groups and reorder particle arrays by cluster.
bool kmeans(ParticleSystem* sys, int* clustersP, int* clusters_size, int k,
            int n_threads);

#endif
