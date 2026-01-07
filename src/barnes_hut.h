#ifndef BARNES_HUT_H
#define BARNES_HUT_H

#include "ds.h"
#include "kmeans.h"

#define EPSILON_O 1e-3
#define CHUNK_SIZE 8

/** Barnes-Hut algorithm
 * Based on the clustered particles, uterlize the space limitation to improve the cache performance.
 * - Build the local tree for each cluster, parallelize the instertion of the particles => no significant improvement.
 * - Or build the global tree
 *  -------------------------------------------------- */
void barnes_hut(double* pos_x, double* pos_y, double* mass, int N, int cluster[][N], double* region, int* clusters_size, int k, double* fx, double* fy, int n_threads, double THETA_MAX);

/** Create a new tree node
 * According the position of the subsquare, determine the region of it.
 * -------------------------------------------------- */
TNode* create_new_TNode(int index, double LB, double RB, double DB, double UB);

/** Insert the particle into the tree
 * If the node represents a particle, then split the node and insert the particle into the corresponding child node.
 * If the node is representing a big object, then update the massive and the centroid of the node.
 * - determine the region the particle belongs to and insert it into the corresponding child node.
 * - if the child node is NULL, then create a new node and insert the particle into it.
 * - else, recursively insert the particle into the child node.
 * -------------------------------------------------- */
int insert(TNode* tNode, double pos_x, double pos_y, double mass, int PID);

/** Compute the force between the particle and the node
 * If the node is a leaf node, directly calculate the force.
 * If the node is a big object, calculate the theta to check if the node is far enough, then calculate the force.
 * - calculate the distance between the particle and the node.
 * - calculate the force factor and update the force.
 * -------------------------------------------------- */
void compute_force(double pos_x, double pos_y, double mass, int PID, TNode* tNode, double* fx, double* fy, int N, double THETA_MAX);

/** Destroy the tree recursively
 *  -------------------------------------------------- */
void destroy(TNode* root);

#endif