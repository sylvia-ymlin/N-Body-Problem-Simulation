#ifndef DS_H
#define DS_H

#include <stddef.h>
#include <stdlib.h>

typedef struct TNode {
    double x_min, x_max, y_min, y_max;
    double pos_x, pos_y, mass;
    int particle_idx;
    struct TNode* child[4];
} TNode;

typedef struct {
    TNode* buffer;
    size_t capacity;
    size_t size;
} NodeArena;

static inline void init_arena(NodeArena* a, size_t cap) {
    a->buffer = (TNode*)malloc(cap * sizeof(TNode));
    a->capacity = cap;
    a->size = 0;
}

static inline TNode* arena_alloc(NodeArena* a) {
    if (!a || a->size >= a->capacity)
        return NULL;
    return &a->buffer[a->size++];
}

static inline void reset_arena(NodeArena* a) {
    a->size = 0;
}
static inline void free_arena(NodeArena* a) {
    if (a && a->buffer) {
        free(a->buffer);
        a->buffer = NULL;
    }
}

#endif
