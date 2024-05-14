#include "kmeans.h"
#include "ds.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ITERATIONS 50

typedef struct {
    double ctr_x;
    double ctr_y;
    int count;
} CNode;

static bool reorder_array(double* arr, const int* order, int N);
static bool reorder_by_clusters(ParticleSystem* sys, const int* clustersP);
static bool converged(CNode* clusters, double* old_clusters_ctr_x,
                      double* old_clusters_ctr_y, int iterations, int k);
static void get_centroids(ParticleSystem* sys, CNode* clusters, int* labels,
                          int k, int N);
static void assign_labels(ParticleSystem* sys, CNode* clusters, int* labels,
                          int k, int N, int n_threads);

bool kmeans(ParticleSystem* sys, int* clustersP, int* clusters_size, int k,
            int n_threads) {
    int N = sys->N;
    CNode clusters[k];
    int* labels = (int*)malloc(N * sizeof(int));
    if (!labels) {
        fprintf(stderr, "Memory allocation failed for kmeans labels.\n");
        exit(1);
    }

    for (int i = 0; i < k; i++) {
        clusters[i].count = 0;
        clusters[i].ctr_x = sys->pos_x[i];
        clusters[i].ctr_y = sys->pos_y[i];
    }

    double old_clusters_ctr_x[k];
    double old_clusters_ctr_y[k];

    int iterations = 0;
    do {
        for (int i = 0; i < k; i++) {
            old_clusters_ctr_x[i] = clusters[i].ctr_x;
            old_clusters_ctr_y[i] = clusters[i].ctr_y;
        }
        iterations++;
        assign_labels(sys, clusters, labels, k, N, n_threads);
        get_centroids(sys, clusters, labels, k, N);
    } while (!converged(clusters, old_clusters_ctr_x, old_clusters_ctr_y,
                        iterations, k));

    int* cnt = (int*)calloc(k, sizeof(int));
    int* offsets = (int*)malloc(k * sizeof(int));

    if (!cnt || !offsets) {
        free(labels);
        if (cnt)
            free(cnt);
        if (offsets)
            free(offsets);
        return false;
    }

    for (int i = 0; i < k; i++)
        clusters_size[i] = 0;
    for (int i = 0; i < N; i++)
        clusters_size[labels[i]]++;

    offsets[0] = 0;
    for (int i = 1; i < k; i++)
        offsets[i] = offsets[i - 1] + clusters_size[i - 1];

    for (int i = 0; i < N; i++) {
        int c = labels[i];
        if (offsets[c] + cnt[c] < N) {
            clustersP[offsets[c] + cnt[c]] = i;
        }
        cnt[c]++;
    }

    if (!reorder_by_clusters(sys, clustersP)) {
        free(cnt);
        free(offsets);
        free(labels);
        return false;
    }

    free(cnt);
    free(offsets);
    free(labels);
    return true;
}

static bool reorder_by_clusters(ParticleSystem* sys, const int* clustersP) {
    int N = sys->N;
    return reorder_array(sys->pos_x, clustersP, N) &&
           reorder_array(sys->pos_y, clustersP, N) &&
           reorder_array(sys->mass, clustersP, N) &&
           reorder_array(sys->vx, clustersP, N) &&
           reorder_array(sys->vy, clustersP, N) &&
           reorder_array(sys->fx, clustersP, N) &&
           reorder_array(sys->fy, clustersP, N);
}

static bool converged(CNode* clusters, double* old_clusters_ctr_x,
                      double* old_clusters_ctr_y, int iterations, int k) {
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

static void get_centroids(ParticleSystem* sys, CNode* clusters, int* labels,
                          int k, int N) {
    for (int i = 0; i < k; i++) {
        clusters[i].ctr_x = 0;
        clusters[i].ctr_y = 0;
        clusters[i].count = 0;
    }

    for (int i = 0; i < N; i++) {
        int label = labels[i];
        clusters[label].ctr_x += sys->pos_x[i];
        clusters[label].ctr_y += sys->pos_y[i];
        clusters[label].count++;
    }

    for (int i = 0; i < k; i++) {
        if (clusters[i].count == 0) {
            if (i < N) {
                clusters[i].ctr_x = sys->pos_x[i];
                clusters[i].ctr_y = sys->pos_y[i];
            }
        } else {
            clusters[i].ctr_x /= clusters[i].count;
            clusters[i].ctr_y /= clusters[i].count;
        }
    }
}

static void assign_labels(ParticleSystem* sys, CNode* clusters, int* labels,
                          int k, int N, int n_threads) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads)
#endif
    for (int i = 0; i < N; i++) {
        double min_dist = INFINITY;
        int label = 0;
        for (int j = 0; j < k; j++) {
            double dx = sys->pos_x[i] - clusters[j].ctr_x;
            double dy = sys->pos_y[i] - clusters[j].ctr_y;
            double dist = dx * dx + dy * dy;
            if (dist < min_dist) {
                min_dist = dist;
                label = j;
            }
        }
        labels[i] = label;
    }
}

static bool reorder_array(double* arr, const int* order, int N) {
    double* tmp = (double*)malloc((size_t)N * sizeof(double));
    if (!tmp) {
        fprintf(stderr, "Memory allocation failed for kmeans reorder buffer.\n");
        return false;
    }

    for (int i = 0; i < N; i++)
        tmp[i] = arr[order[i]];

    memcpy(arr, tmp, (size_t)N * sizeof(double));
    free(tmp);
    return true;
}
