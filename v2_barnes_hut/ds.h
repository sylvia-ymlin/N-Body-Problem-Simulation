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

#endif
