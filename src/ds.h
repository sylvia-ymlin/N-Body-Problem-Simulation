#ifndef DS_H
#define DS_H

#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


/* Maximum iterations for k-means clustering. */
#define MAX_ITERATIONS 50

/* Structure to represent a tree node
 * - bounding box: LB, RB, DB, UB
 * - child nodes: child[4]
 * - particle/object properties: pos_x, pos_y, mass, index
 * - PID: index of particles or big object (PID = -1)
 */
typedef struct TNode {
  double LB, RB, DB, UB;
  struct TNode *child[4];
  double pos_x;
  double pos_y;
  double mass;
  int PID;
} TNode;

/* Structure to represent a cluster node. */
typedef struct CNode {
  double ctr_x;
  double ctr_y;
  int count; /* Number of particles in the cluster*/
} CNode;

/* Memory Arena for fast TNode allocation */
typedef struct NodeArena {
  TNode *nodes;    // Pre-allocated buffer
  size_t capacity; // Total number of nodes
  size_t used;     // Current index
} NodeArena;

#endif