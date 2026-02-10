#include "kmeans.h"

// K-Means Load Balancing: groups particles into spatially coherent clusters
// to maximize cache reuse and balance OpenMP workloads.
bool kmeans(double *pos_x, double *pos_y, int N, int *clustersP,
            int *clusters_size, int k, int n_threads) {

  CNode clusters[k];
  int *labels = (int *)malloc(N * sizeof(int));
  if (!labels) {
    fprintf(stderr, "Memory allocation failed for kmeans labels.\n");
    exit(1);
  }

  double old_clusters_ctr_x[k];
  double old_clusters_ctr_y[k];

  /* Initialize the clusters centroids. */
  for (int i = 0; i < k; i++) {
    clusters[i].count = 0;
    clusters[i].ctr_x = pos_x[i];
    clusters[i].ctr_y = pos_y[i];
  }

  /* Iterate until the clusters centroids converge */
  int iterations = 0;
  while (!converged(clusters, old_clusters_ctr_x, old_clusters_ctr_y,
                    iterations, k)) {
    /* Save the old clusters centroids. */
    for (int i = 0; i < k; i++) {
      old_clusters_ctr_x[i] = clusters[i].ctr_x;
      old_clusters_ctr_y[i] = clusters[i].ctr_y;
    }

    iterations++;

    /* Assign the particles to the clusters. */
    assignLabels(pos_x, pos_y, clusters, labels, k, N, n_threads);

    /* Update clusters */
    getCentroids(pos_x, pos_y, clusters, labels, k, N);
  }

  /* Group particles into clusters (CSR-like packing) */
  int *cnt = (int *)calloc(k, sizeof(int));
  int *offsets = (int *)malloc(k * sizeof(int));

  // 1. Count sizes
  for (int i = 0; i < k; i++)
    clusters_size[i] = 0;
  for (int i = 0; i < N; i++)
    clusters_size[labels[i]]++;

  // 2. Compute offsets
  offsets[0] = 0;
  for (int i = 1; i < k; i++)
    offsets[i] = offsets[i - 1] + clusters_size[i - 1];

  // 3. Fill packed array
  for (int i = 0; i < N; i++) {
    int c = labels[i];
    clustersP[offsets[c] + cnt[c]] = i;
    cnt[c]++;
  }

  free(cnt);
  free(offsets);
  free(labels);
  return true;
}

bool converged(CNode *clusters, double *old_clusters_ctr_x,
               double *old_clusters_ctr_y, int iterations, int k) {
  if (iterations > MAX_ITERATIONS) {
    return true;
  }

  for (int i = 0; i < k; i++) {
    if (fabs(clusters[i].ctr_x - old_clusters_ctr_x[i]) > 1e-5 ||
        fabs(clusters[i].ctr_y - old_clusters_ctr_y[i]) > 1e-5) {
      return false;
    }
  }

  return true;
}

void getCentroids(double *pos_x, double *pos_y, CNode *clusters, int *labels,
                  int k, int N) {
  /* reset the clusters. */
  for (int i = 0; i < k; i++) {
    clusters[i].ctr_x = 0;
    clusters[i].ctr_y = 0;
    clusters[i].count = 0;
  }

  /* update the clusters. */
  for (int i = 0; i < N; i++) {
    clusters[labels[i]].ctr_x += pos_x[i];
    clusters[labels[i]].ctr_y += pos_y[i];
    clusters[labels[i]].count++;
  }

  /* calculate the centroids. */
  for (int i = 0; i < k; i++) {
    if (clusters[i].count == 0) {
      clusters[i].ctr_x = pos_x[i];
      clusters[i].ctr_y = pos_y[i];
    } else {
      clusters[i].ctr_x /= clusters[i].count;
      clusters[i].ctr_y /= clusters[i].count;
    }
  }
}

void assignLabels(double *pos_x, double *pos_y, CNode *clusters, int *labels,
                  int k, int N, int n_threads) {
#pragma omp parallel for num_threads(n_threads)
  for (int i = 0; i < N; i++) {
    double min_dist = INFINITY;
    int label = 0;
    for (int j = 0; j < k; j++) {
      double dist =
          (pos_x[i] - clusters[j].ctr_x) * (pos_x[i] - clusters[j].ctr_x) +
          (pos_y[i] - clusters[j].ctr_y) * (pos_y[i] - clusters[j].ctr_y);
      if (dist < min_dist) {
        min_dist = dist;
        label = j;
      }
    }
    labels[i] = label;
  }
}
