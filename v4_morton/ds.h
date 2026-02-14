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

typedef struct {
  TNode *buffer;
  size_t size;
  size_t used;
} NodeArena;

static inline void init_arena(NodeArena *arena, size_t size) {
  arena->buffer = (TNode *)malloc(size * sizeof(TNode));
  arena->size = size;
  arena->used = 0;
}

static inline void free_arena(NodeArena *arena) {
  if (arena->buffer) {
    free(arena->buffer);
    arena->buffer = NULL;
  }
}

static inline void reset_arena(NodeArena *arena) {
  arena->used = 0;
}

static inline TNode *arena_alloc(NodeArena *arena) {
  if (arena && arena->used < arena->size) {
    return &arena->buffer[arena->used++];
  }
  return NULL; // Fallback or error
}

#endif
