#ifndef DS_H
#define DS_H

#include <stdlib.h>

typedef struct TNode {
  double LB, RB, DB, UB;
  struct TNode *child[4];
  double pos_x;
  double pos_y;
  double mass;
  int PID;
} TNode;

typedef struct NodeArena {
  TNode *nodes;
  size_t capacity;
  size_t used;
} NodeArena;

static inline void init_arena(NodeArena *arena, size_t capacity) {
  arena->nodes = (TNode *)malloc(capacity * sizeof(TNode));
  arena->capacity = capacity;
  arena->used = 0;
}

static inline void free_arena(NodeArena *arena) {
  free(arena->nodes);
}

static inline void reset_arena(NodeArena *arena) {
  arena->used = 0;
}

static inline TNode *arena_alloc(NodeArena *arena) {
  if (arena->used < arena->capacity) {
    return &arena->nodes[arena->used++];
  }
  return NULL;
}

#endif
