#include "../utils/ds.h"
#include "kernels.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static NodeArena arena = {NULL, 0, 0};
static int arena_N = 0;

static int is_leaf(TNode *node) { return node->child[0] == NULL; }

static TNode *create_node(double LB, double RB, double DB, double UB,
                          NodeArena *a) {
  TNode *node = arena_alloc(a);
  if (!node) {
    fprintf(stderr, "Error: Arena out of memory!\n");
    exit(1);
  }
  node->LB = LB;
  node->RB = RB;
  node->DB = DB;
  node->UB = UB;
  node->pos_x = 0;
  node->pos_y = 0;
  node->mass = 0;
  node->PID = -1;
  for (int i = 0; i < 4; i++)
    node->child[i] = NULL;
  return node;
}

static void insert(TNode *node, int idx, ParticleSystem *sys, NodeArena *a) {
  double px = sys->pos_x[idx];
  double py = sys->pos_y[idx];
  double mass = sys->mass[idx];

  if (node->mass == 0 && is_leaf(node) && node->PID == -1) {
    node->mass = mass;
    node->pos_x = px;
    node->pos_y = py;
    node->PID = idx;
    return;
  }

  if (is_leaf(node)) {
    if (node->PID != -1) {
      double dx = fabs(px - sys->pos_x[node->PID]);
      double dy = fabs(py - sys->pos_y[node->PID]);
      if (dx < 1e-9 && dy < 1e-9) {
        double total_mass = node->mass + mass;
        node->pos_x = (node->pos_x * node->mass + px * mass) / total_mass;
        node->pos_y = (node->pos_y * node->mass + py * mass) / total_mass;
        node->mass = total_mass;
        return;
      }

      int old_idx = node->PID;
      node->PID = -1;
      insert(node, old_idx, sys, a);
    }
  }

  double mid_x = (node->LB + node->RB) * 0.5;
  double mid_y = (node->DB + node->UB) * 0.5;
  int quad = 0;
  if (px > mid_x)
    quad += 1;
  if (py > mid_y)
    quad += 2;

  if (!node->child[quad]) {
    double nx_min = (quad & 1) ? mid_x : node->LB;
    double nx_max = (quad & 1) ? node->RB : mid_x;
    double ny_min = (quad & 2) ? mid_y : node->DB;
    double ny_max = (quad & 2) ? node->UB : mid_y;
    node->child[quad] = create_node(nx_min, nx_max, ny_min, ny_max, a);
  }
  insert(node->child[quad], idx, sys, a);

  double total_mass = node->mass + mass;
  node->pos_x = (node->pos_x * node->mass + px * mass) / total_mass;
  node->pos_y = (node->pos_y * node->mass + py * mass) / total_mass;
  node->mass = total_mass;
}

static void compute_force_rec(TNode *node, int idx, ParticleSystem *sys,
                              double *fx, double *fy, double THETA_MAX,
                              double G) {
  if (!node || node->mass == 0 || node->PID == idx)
    return;

  double dx = node->pos_x - sys->pos_x[idx];
  double dy = node->pos_y - sys->pos_y[idx];
  double r = sqrt(dx * dx + dy * dy + 1e-6);
  double s = node->RB - node->LB;

  if (is_leaf(node) || (s / r < THETA_MAX)) {
    double f = G * sys->mass[idx] * node->mass / (r * r * r);
    *fx += f * dx;
    *fy += f * dy;
  } else {
    for (int i = 0; i < 4; i++) {
      compute_force_rec(node->child[i], idx, sys, fx, fy, THETA_MAX, G);
    }
  }
}

void compute_force_v3_arena(ParticleSystem *sys, KernelConfig *config) {
  if (arena.buffer == NULL || arena_N != sys->N) {
    if (arena.buffer)
      free_arena(&arena);
    // Estimate nodes needed. 2.0 * N is usually sufficient for quadtree, but
    // for safety 4.0 or more. v3 used 1000 * N, which is huge overkill but
    // safe. I'll use 100 * N to be safe but not insane.
    init_arena(&arena, sys->N * 100);
    arena_N = sys->N;
  }
  reset_arena(&arena);

  int N = sys->N;
  double x_min = sys->pos_x[0], x_max = sys->pos_x[0];
  double y_min = sys->pos_y[0], y_max = sys->pos_y[0];

  for (int i = 1; i < N; i++) {
    if (sys->pos_x[i] < x_min)
      x_min = sys->pos_x[i];
    if (sys->pos_x[i] > x_max)
      x_max = sys->pos_x[i];
    if (sys->pos_y[i] < y_min)
      y_min = sys->pos_y[i];
    if (sys->pos_y[i] > y_max)
      y_max = sys->pos_y[i];
  }

  TNode *root = create_node(x_min, x_max, y_min, y_max, &arena);
  for (int i = 0; i < N; i++) {
    insert(root, i, sys, &arena);
  }

  double G = 100.0 / N;
  double theta = config->theta_max;

  for (int i = 0; i < N; i++) {
    sys->fx[i] = 0;
    sys->fy[i] = 0;
    compute_force_rec(root, i, sys, &sys->fx[i], &sys->fy[i], theta, G);
  }
}
