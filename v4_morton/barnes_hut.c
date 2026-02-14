#include "barnes_hut.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static TNode *create_node(double LB, double RB, double DB, double UB, NodeArena *arena) {
  TNode *node;
  if (arena) {
    node = arena_alloc(arena);
  } else {
    node = (TNode *)malloc(sizeof(TNode));
  }
  
  if (!node) {
      fprintf(stderr, "Memory allocation failed!\n");
      exit(1);
  }
  node->LB = LB;
  node->RB = RB;
  node->DB = DB;
  node->UB = UB;
  node->mass = 0;
  node->pos_x = 0;
  node->pos_y = 0;
  node->PID = -1;
  for (int i = 0; i < 4; i++)
    node->child[i] = NULL;
  return node;
}

static void free_tree(TNode *node) {
  if (!node) return;
  for (int i = 0; i < 4; i++) {
    if (node->child[i]) {
      free_tree(node->child[i]);
    }
  }
  free(node);
}

static void insert(TNode *node, int idx, double *px, double *py, double *mass, NodeArena *arena) {
  if (node->PID == -1 && node->mass == 0) {
    node->PID = idx;
    node->pos_x = px[idx];
    node->pos_y = py[idx];
    node->mass = mass[idx];
    return;
  }

  if (node->PID != -1) {
    double dx = fabs(px[idx] - px[node->PID]);
    double dy = fabs(py[idx] - py[node->PID]);
    if (dx < 1e-9 && dy < 1e-9) {
       double total_mass = node->mass + mass[idx];
       node->pos_x = (node->pos_x * node->mass + px[idx] * mass[idx]) / total_mass;
       node->pos_y = (node->pos_y * node->mass + py[idx] * mass[idx]) / total_mass;
       node->mass = total_mass;
       return;
    }
    int old_idx = node->PID;
    node->PID = -1;
    insert(node, old_idx, px, py, mass, arena);
  }

  double mid_x = (node->LB + node->RB) * 0.5;
  double mid_y = (node->DB + node->UB) * 0.5;
  int quad = (px[idx] > mid_x) + ((py[idx] > mid_y) << 1);

  if (!node->child[quad]) {
    double nLB = (quad & 1) ? mid_x : node->LB;
    double nRB = (quad & 1) ? node->RB : mid_x;
    double nDB = (quad & 2) ? mid_y : node->DB;
    double nUB = (quad & 2) ? node->UB : mid_y;
    node->child[quad] = create_node(nLB, nRB, nDB, nUB, arena);
  }
  insert(node->child[quad], idx, px, py, mass, arena);

  double total_mass = node->mass + mass[idx];
  node->pos_x = (node->pos_x * node->mass + px[idx] * mass[idx]) / total_mass;
  node->pos_y = (node->pos_y * node->mass + py[idx] * mass[idx]) / total_mass;
  node->mass = total_mass;
}

static void compute_force(TNode *node, int idx, double *px, double *py,
                          double *mass, double *fx, double *fy,
                          double THETA_MAX, double G) {
  if (!node || node->mass == 0 || (node->PID == idx && node->PID != -1))
    return;

  double dx = node->pos_x - px[idx];
  double dy = node->pos_y - py[idx];
  double r = sqrt(dx * dx + dy * dy + 1e-6);
  double s = node->RB - node->LB;

  if (node->PID != -1 || (s / r < THETA_MAX)) {
    double f = G * mass[idx] * node->mass / (r * r * r);
    *fx += f * dx;
    *fy += f * dy;
  } else {
    for (int i = 0; i < 4; i++) {
      if (node->child[i])
        compute_force(node->child[i], idx, px, py, mass, fx, fy, THETA_MAX, G);
    }
  }
}

void barnes_hut(double *px, double *py, double *mass, int N, double *fx,
                double *fy, double THETA_MAX, int use_arena) {
  
  for (int i = 0; i < N; i++) {
    fx[i] = 0;
    fy[i] = 0;
  }

  double x_min = px[0], x_max = px[0], y_min = py[0], y_max = py[0];
  for (int i = 1; i < N; i++) {
    if (px[i] < x_min)
      x_min = px[i];
    if (px[i] > x_max)
      x_max = px[i];
    if (py[i] < y_min)
      y_min = py[i];
    if (py[i] > y_max)
      y_max = py[i];
  }
  
  // Square region
  double dx = x_max - x_min;
  double dy = y_max - y_min;
  double d = (dx > dy) ? dx : dy;
  x_max = x_min + d;
  y_max = y_min + d;
  
  // Padding
  x_min -= d * 0.05;
  x_max += d * 0.05;
  y_min -= d * 0.05;
  y_max += d * 0.05;

  NodeArena arena_struct;
  NodeArena *arena = NULL;
  
  if (use_arena) {
    // Heuristic: 2N nodes is usually enough for Barnes-Hut
    // But let's be safe with 4N or even more to avoid realloc (we don't have realloc in simple arena)
    init_arena(&arena_struct, N * 10); // Safe upper bound
    arena = &arena_struct;
  }

  TNode *root = create_node(x_min, x_max, y_min, y_max, arena);
  
  for (int i = 0; i < N; i++) {
    insert(root, i, px, py, mass, arena);
  }

  double G = 100.0 / N;
  for (int i = 0; i < N; i++) {
    compute_force(root, i, px, py, mass, &fx[i], &fy[i], THETA_MAX, G);
  }

  // IMPORTANT: Free the tree
  if (use_arena) {
    free_arena(arena);
  } else {
    free_tree(root);
  }
}
