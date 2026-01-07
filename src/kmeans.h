#ifndef KMEANS_H
#define KMEANS_H

#include "ds.h"

/* Function: kmeans
 * --------------------------------------------------
 * This function implements the K-Means algorithm, which is a clustering algorithm.
 * It clusters the particles into k clusters based on the distance to the cluster centroids.
 * The function returns the groups of particles in the clusters and the size of each cluster.
 */
bool kmeans(double* pos_x, double* pos_y, int N, int cluster[][N], int* clusters_size, int k, int n_threads);

/* Function: converged
 * --------------------------------------------------
 * This function checks if the clusters centroids have converged.
 * K Means terminates when the clusters centroids do not change or achieve a maximum number of iterations. 
 */
bool converged(CNode* clusters, double* old_clusters_pos_x, double* old_clusters_pos_y, int iterations, int k);

/* Function: get the centroids
 * --------------------------------------------------
 * This function initializes the clusters centroids when the corresponding clusters are empty.
 * Or calculates the centroids as the geometric center of the particles in the cluster. 
 */
void getCentroids(double* pos_x, double* pos_y, CNode* clusters, int* labels, int k, int N);

/* Function: assign labels for the particles
 * --------------------------------------------------
 * This function assigns the particles to the clusters based on the distance to the cluster centroids.
 */
void assignLabels(double* pos_x, double* pos_y, CNode* clusters, int* labels, int k, int N, int n_threads);

#endif