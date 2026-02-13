#ifndef DS_H
#define DS_H

#include <stdlib.h>

typedef struct TNode {
  double mass;
  double x, y;
  double x_min, x_max, y_min, y_max;
  struct TNode *children[4];
  int is_leaf;
  int particle_idx;
} TNode;

typedef struct {
  TNode *buffer;
  int size;
  int used;
} NodeArena;

static inline void init_arena(NodeArena *arena, int size) {
  arena->buffer = (TNode *)malloc(size * sizeof(TNode));
  arena->size = size;
  arena->used = 0;
}

static inline void free_arena(NodeArena *arena) {
  free(arena->buffer);
}

static inline void reset_arena(NodeArena *arena) {
  arena->used = 0;
}

static inline TNode *arena_alloc(NodeArena *arena) {
  if (arena->used < arena->size) {
    return &arena->buffer[arena->used++];
  }
  return NULL;
}

#endif
