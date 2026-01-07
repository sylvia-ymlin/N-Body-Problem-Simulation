#include "kmeans.h"

/* Function: kmeans
 * --------------------------------------------------
 * This function implements the K-Means algorithm, which is a clustering algorithm.
 * It clusters the particles into k clusters based on the distance to the cluster centroids.
 * The function returns the groups of particles in the clusters and the size of each cluster through pointers.
 */
bool kmeans(double* pos_x, double* pos_y, int N, int clustersP[][N], int* clusters_size, int k, int n_threads) {

    CNode clusters[k];
    int labels[N];
    double old_clusters_ctr_x[k];
    double old_clusters_ctr_y[k];

    /* Initialize the clusters centroids. */
    for (int i = 0; i < k; i++) {
        clusters[i].count = 0;
        clusters[i].ctr_x = pos_x[i];
        clusters[i].ctr_y = pos_y[i];
    }

    /* Iterate until the clusters centroids converge */
    int interations = 0;
    while (!converged(clusters, old_clusters_ctr_x, old_clusters_ctr_y, interations, k)) {
        /* Save the old clusters centroids. */
        for (int i = 0; i < k; i++) {
            old_clusters_ctr_x[i] = clusters[i].ctr_x;
            old_clusters_ctr_y[i] = clusters[i].ctr_y;
        }

        interations++;

        /* Assign the particles to the clusters. */
        assignLabels(pos_x, pos_y, clusters, labels, k, N, n_threads);

        /* Update the clusters centroids. */
        getCentroids(pos_x, pos_y, clusters, labels, k, N);
    }

    /* Group the particles into the clusters. */
    for (int i = 0; i < k; i++) {
        clusters_size[i] = 0;
    }
    for (int i = 0; i < N; i++) {
        clustersP[labels[i]][clusters_size[labels[i]]++] = i;
    }
    return true;
}

/* Function: converged
 * --------------------------------------------------
 * This function checks if the clusters centroids have converged.
 * K Means terminates when the clusters centroids do not change or achieve a maximum number of iterations. 
 */
bool converged(CNode* clusters, double* old_clusters_ctr_x, double* old_clusters_ctr_y, int iterations, int k) {
    if (iterations > MAX_ITERATIONS) {
        return true;
    }

    for (int i = 0; i < k; i++) {
        if (clusters[i].ctr_x != old_clusters_ctr_x[i] || clusters[i].ctr_y != old_clusters_ctr_y[i]) {
            return false;
        }
    }

    return true;
}

/* Function: get the centroids
 * --------------------------------------------------
 * This function initializes the clusters centroids when the corresponding clusters are empty.
 * Or calculates the centroids as the geometric center of the particles in the cluster. 
 */
void getCentroids(double* pos_x, double* pos_y, CNode* clusters, int* labels, int k, int N) {
    /* reset the clusters. */
    for (int i = 0; i < k; i++) {
        clusters[i].ctr_x = 0;
        clusters[i].ctr_y = 0;
        clusters[i].count = 0;
    }

    /* sum the positions of the particles in the same cluster and count the number of particles in the cluster. */
    for (int i = 0; i < N; i++) {
        clusters[labels[i]].ctr_x += pos_x[i];
        clusters[labels[i]].ctr_y += pos_y[i];
        clusters[labels[i]].count++;
    }

    /* calculate the new centroids or initialize the centroids if the cluster is empty. */
    for (int i = 0; i < k; i++) {
        clusters[i].ctr_x = clusters[i].count == 0 ? pos_x[i] : clusters[i].ctr_x / clusters[i].count;
        clusters[i].ctr_y = clusters[i].count == 0 ? pos_y[i] : clusters[i].ctr_y / clusters[i].count;
    }
}

/* Function: assign labels for the particles
 * --------------------------------------------------
 * This function assigns the particles to the clusters based on the distance to the cluster centroids.
 */
void assignLabels(double* pos_x, double* pos_y, CNode* clusters, int* labels, int k, int N, int n_threads) {
    /* Paralize the computation of loss function. */
    # ifdef _OPENMP
        #pragma omp parallel for schedule(static, 6)  num_threads(n_threads)
    # endif
    for (int i = 0; i < N; i++) {
        double loss = INFINITY, cur_loss;
        for (int j = 0; j < k; j++) {
            cur_loss = pow(pos_x[i] - clusters[j].ctr_x, 2) +
                       pow(pos_y[i] - clusters[j].ctr_y, 2);
            if (cur_loss < loss) {
                loss = cur_loss;
                labels[i] = j;
            }
        }
    }
}