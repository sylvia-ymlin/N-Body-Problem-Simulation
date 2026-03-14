#include "simulation.h"
#include "src/ds.h"
#include "src/morton.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * v4: Morton Ordering - Spatial locality
 *
 * Maps 2D spatial proximity to 1D memory contiguity via Z-order curves.
 * Forces CPU cache alignment for recursive tree traversals.
 */

static NodeArena arena = {NULL, 0, 0};
static int arena_N = 0;

static inline int is_leaf(TNode *node) {
  return node->child[0] == NULL && node->child[1] == NULL &&
         node->child[2] == NULL && node->child[3] == NULL;
}

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
  double px = sys->pos_x[idx], py = sys->pos_y[idx], m = sys->mass[idx];

  if (node->PID == -1 && is_leaf(node) && node->mass == 0) {
    node->PID = idx;
    node->mass = m;
    node->pos_x = px;
    node->pos_y = py;
    return;
  }

  if (is_leaf(node)) {
    int old = node->PID;
    if (old != -1) {
      if (fabs(px - sys->pos_x[old]) < 1e-9 &&
          fabs(py - sys->pos_y[old]) < 1e-9) {
        double mt = node->mass + m;
        node->pos_x = (node->pos_x * node->mass + px * m) / mt;
        node->pos_y = (node->pos_y * node->mass + py * m) / mt;
        node->mass = mt;
        return;
      }
      node->PID = -1;
      double old_px = node->pos_x, old_py = node->pos_y, old_m = node->mass;
      node->mass = 0;
      node->pos_x = 0;
      node->pos_y = 0;

      double mx = (node->LB + node->RB) * 0.5, my = (node->DB + node->UB) * 0.5;
      int oq = (old_px > mx) + 2 * (old_py > my);
      node->child[oq] =
          create_node((oq & 1) ? mx : node->LB, (oq & 1) ? node->RB : mx,
                      (oq & 2) ? my : node->DB, (oq & 2) ? node->UB : my, a);
      insert(node->child[oq], old, sys, a);
    }
  }

  double mx = (node->LB + node->RB) * 0.5, my = (node->DB + node->UB) * 0.5;
  int q = (px > mx) + 2 * (py > my);
  if (!node->child[q]) {
    node->child[q] =
        create_node((q & 1) ? mx : node->LB, (q & 1) ? node->RB : mx,
                    (q & 2) ? my : node->DB, (q & 2) ? node->UB : my, a);
  }
  insert(node->child[q], idx, sys, a);

  double mt = node->mass + m;
  node->pos_x = (node->pos_x * node->mass + px * m) / mt;
  node->pos_y = (node->pos_y * node->mass + py * m) / mt;
  node->mass = mt;
}

static void compute_force_rec(TNode *node, int idx, ParticleSystem *sys,
                              double *fx, double *fy, double THETA, double G) {
  if (!node || node->mass == 0 || node->PID == idx)
    return;

  double dx = node->pos_x - sys->pos_x[idx];
  double dy = node->pos_y - sys->pos_y[idx];
  double r = sqrt(dx * dx + dy * dy + 1e-6);
  double s = node->RB - node->LB;

  if (is_leaf(node) || (s / r < THETA)) {
    double f = G * sys->mass[idx] * node->mass / (r * r * r);
    *fx += f * dx;
    *fy += f * dy;
  } else {
    for (int i = 0; i < 4; i++)
      compute_force_rec(node->child[i], idx, sys, fx, fy, THETA, G);
  }
}

void compute_force_v4_morton(ParticleSystem *sys, KernelConfig *config) {
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

  // Ensure square bounding box and add padding for sorting
  double dx = x_max - x_min;
  double dy = y_max - y_min;
  double d = (dx > dy) ? dx : dy;
  x_max = x_min + d;
  y_max = y_min + d;
  x_min -= d * 0.05;
  x_max += d * 0.05;
  y_min -= d * 0.05;
  y_max += d * 0.05;

  // STEP 1: Sort by Morton Code (Z-order)
  z_order_sort(sys, x_min, x_max, y_min, y_max);

  // STEP 2: Build Tree using Arena
  if (arena.buffer == NULL || arena_N != sys->N) {
    if (arena.buffer)
      free_arena(&arena);
    init_arena(&arena, sys->N * 100);
    arena_N = sys->N;
  }
  reset_arena(&arena);

  TNode *root = create_node(x_min, x_max, y_min, y_max, &arena);
  for (int i = 0; i < N; i++)
    insert(root, i, sys, &arena);

  // STEP 3: Compute Forces
  double G = 100.0 / N;
  double theta = config->theta_max;

  for (int i = 0; i < N; i++) {
    sys->fx[i] = 0;
    sys->fy[i] = 0;
    compute_force_rec(root, i, sys, &sys->fx[i], &sys->fy[i], theta, G);
  }
}
