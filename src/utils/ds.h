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
  int PID; // Particle ID, or -1 for internal node
} TNode;

/* Structure to represent a cluster node. */
typedef struct CNode {
  double ctr_x;
  double ctr_y;
  int count; /* Number of particles in the cluster*/
} CNode;

typedef struct {
  TNode *buffer;
  size_t size;
  size_t used;
} NodeArena;

static inline void init_arena(NodeArena *arena, size_t size) {
  arena->buffer = (TNode *)malloc(size * sizeof(TNode));
  if (arena->buffer == NULL) {
    fprintf(stderr, "Error: Failed to allocate arena of size %zu\n", size);
    exit(1);
  }
  arena->size = size;
  arena->used = 0;
}

static inline void free_arena(NodeArena *arena) {
  if (arena->buffer) {
    free(arena->buffer);
    arena->buffer = NULL;
  }
}

static inline void reset_arena(NodeArena *arena) { arena->used = 0; }

static inline TNode *arena_alloc(NodeArena *arena) {
  if (arena && arena->used < arena->size) {
    return &arena->buffer[arena->used++];
  }
  return NULL; // Fallback or error
}

#endif
