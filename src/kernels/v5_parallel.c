#include "../utils/ds.h"
#include "../utils/kmeans.h"
#include "../utils/morton.h"
#include "kernels.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

static int *clusters = NULL;
static int *clusters_size = NULL;
static int allocated_N = 0;
static int allocated_k = 0;
static int call_count = 0;

static NodeArena arena = {NULL, 0, 0};
static int arena_N = 0;

static TNode *create_new_TNode(int index, double LB, double RB, double DB,
                               double UB, NodeArena *arena) {
  TNode *new_TNode;
  if (arena) {
    new_TNode = arena_alloc(arena);
    if (!new_TNode) {
      fprintf(stderr, "Error: Arena out of memory inside v5!\n");
      exit(1);
    }
  } else {
    new_TNode = (TNode *)malloc(sizeof(TNode));
  }

  for (int i = 0; i < 4; i++)
    new_TNode->child[i] = NULL;
  new_TNode->PID = -1;
  new_TNode->mass = 0;
  new_TNode->pos_x = 0;
  new_TNode->pos_y = 0;

  double mid_x = 0.5 * (LB + RB);
  double mid_y = 0.5 * (DB + UB);

  if (index == -1) {
    new_TNode->LB = LB;
    new_TNode->RB = RB;
    new_TNode->DB = DB;
    new_TNode->UB = UB;
  } else {
    new_TNode->LB = (index == 2 || index == 3) ? mid_x : LB;
    new_TNode->RB = (index == 2 || index == 3) ? RB : mid_x;
    new_TNode->DB = (index == 1 || index == 3) ? mid_y : DB;
    new_TNode->UB = (index == 1 || index == 3) ? UB : mid_y;
  }
  return new_TNode;
}

static void free_tree_rec(TNode *node) {
  if (!node)
    return;
  for (int i = 0; i < 4; i++)
    free_tree_rec(node->child[i]);
  free(node);
}

// Insert logic from v5
static void insert(TNode *tNode, double pos_x, double pos_y, double mass,
                   int PID, NodeArena *arena) {
  double mid_x = 0.5 * (tNode->LB + tNode->RB);
  double mid_y = 0.5 * (tNode->DB + tNode->UB);
  double width = tNode->RB - tNode->LB;

  if (tNode->PID != -1) {
    if (width < 1e-12 || ((pos_x == tNode->pos_x) && (pos_y == tNode->pos_y))) {
      double new_mass = tNode->mass + mass;
      tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
      tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
      tNode->mass = new_mass;
      return;
    }

    int old_index = (tNode->pos_y > mid_y) + 2 * (tNode->pos_x > mid_x);
    tNode->child[old_index] = create_new_TNode(old_index, tNode->LB, tNode->RB,
                                               tNode->DB, tNode->UB, arena);
    tNode->child[old_index]->PID = tNode->PID;
    tNode->child[old_index]->pos_x = tNode->pos_x;
    tNode->child[old_index]->pos_y = tNode->pos_y;
    tNode->child[old_index]->mass = tNode->mass;

    tNode->PID = -1;
  }

  double new_mass = tNode->mass + mass;
  if (new_mass > 0) {
    tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
    tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
  } else {
    tNode->pos_x = pos_x;
    tNode->pos_y = pos_y;
  }
  tNode->mass = new_mass;

  int index = (pos_y > mid_y) + 2 * (pos_x > mid_x);
  if (tNode->child[index] == NULL) {
    tNode->child[index] = create_new_TNode(index, tNode->LB, tNode->RB,
                                           tNode->DB, tNode->UB, arena);
    tNode->child[index]->pos_x = pos_x;
    tNode->child[index]->pos_y = pos_y;
    tNode->child[index]->mass = mass;
    tNode->child[index]->PID = PID;
  } else {
    insert(tNode->child[index], pos_x, pos_y, mass, PID, arena);
  }
}

static inline void _compute_force(double pos_x, double pos_y, double mass,
                                  TNode *tNode, double *fx, double *fy,
                                  double G) {
  double r_x = tNode->pos_x - pos_x;
  double r_y = tNode->pos_y - pos_y;
  double r_sq = r_x * r_x + r_y * r_y + 1e-12; // Epsilon
  double r_inv = 1.0 / sqrt(r_sq);
  double r_inv3 = r_inv * r_inv * r_inv;
  double force_factor = G * mass * tNode->mass * r_inv3;
  *fx += force_factor * r_x;
  *fy += force_factor * r_y;
}

static void compute_force_single(double pos_x, double pos_y, double mass,
                                 int PID, TNode *root, double *fx, double *fy,
                                 double G, double THETA_MAX) {
  TNode *stack[2048];
  int sp = 0;
  if (root)
    stack[sp++] = root;

  while (sp > 0) {
    TNode *tNode = stack[--sp];
    if (tNode->PID == PID)
      continue;

    double width = (tNode->RB - tNode->LB);
    double dx = pos_x - tNode->pos_x;
    double dy = pos_y - tNode->pos_y;
    double dist_sq = dx * dx + dy * dy;

    if (tNode->PID != -1) {
      _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
    } else {
      if (width * width < THETA_MAX * THETA_MAX * dist_sq) {
        _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
      } else {
        for (int i = 0; i < 4; i++) {
          if (tNode->child[i]) {
            if (sp < 2048)
              stack[sp++] = tNode->child[i];
          }
        }
      }
    }
  }
}

void compute_force_v5_parallel(ParticleSystem *sys, KernelConfig *config) {
  int N = sys->N;
  int k = config->k_clusters;
  // Note: if k <= 0, we run in ablation mode

  // Calculate Bounds (needed for Sort AND Tree)
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

  // Make square and pad (consistent with v4/v5)
  double dx = x_max - x_min;
  double dy = y_max - y_min;
  double d = (dx > dy) ? dx : dy;
  x_max = x_min + d;
  y_max = y_min + d;
  x_min -= d * 0.05;
  x_max += d * 0.05;
  y_min -= d * 0.05;
  y_max += d * 0.05;

  // ABLATION: If k=0, Sort Particles (MORTON)
  if (k <= 0) {
    z_order_sort(sys, x_min, x_max, y_min, y_max);
  } else {
    // 1. Manage Clusters Memory (Only if k > 0)
    if (!clusters || allocated_N < N) {
      if (clusters)
        free(clusters);
      clusters = (int *)malloc(N * sizeof(int));
      allocated_N = N;
    }
    if (!clusters_size || allocated_k < k) {
      if (clusters_size)
        free(clusters_size);
      clusters_size = (int *)malloc(k * sizeof(int));
      allocated_k = k;
    }

    // 2. Run KMeans periodically
    if (call_count == 0 || (call_count % 10 == 0)) {
      kmeans(sys, clusters, clusters_size, k, config->n_threads);
    }
  }
  call_count++;

  // 3. Reset Forces
  for (int i = 0; i < N; i++) {
    sys->fx[i] = 0;
    sys->fy[i] = 0;
  }

  // 4. Build Tree
  NodeArena *a_ptr = NULL;
  if (config->use_arena) {
    if (arena.buffer == NULL || arena_N != N) {
      if (arena.buffer)
        free_arena(&arena);
      init_arena(&arena, N * 10);
      arena_N = N;
    }
    reset_arena(&arena);
    a_ptr = &arena;
  }

  TNode *root = create_new_TNode(-1, x_min, x_max, y_min, y_max, a_ptr);
  root->PID = 0;
  root->pos_x = sys->pos_x[0];
  root->pos_y = sys->pos_y[0];
  root->mass = sys->mass[0];

  for (int i = 1; i < N; i++) {
    insert(root, sys->pos_x[i], sys->pos_y[i], sys->mass[i], i,
           a_ptr); // i is correct PID because we use position in array
  }

  // 5. Compute Forces (Parallel)
  double G = 100.0 / N;

  if (k > 0) {
    // Original K-Means based scheduling
    int *offsets = (int *)malloc(k * sizeof(int));
    offsets[0] = 0;
    for (int i = 1; i < k; i++)
      offsets[i] = offsets[i - 1] + clusters_size[i - 1];

#pragma omp parallel for schedule(dynamic, 8) num_threads(config->n_threads)
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < clusters_size[i]; j++) {
        int PID = clusters[offsets[i] + j];
        compute_force_single(sys->pos_x[PID], sys->pos_y[PID], sys->mass[PID],
                             PID, root, &sys->fx[PID], &sys->fy[PID], G,
                             config->theta_max);
      }
    }
    free(offsets);
  } else {
// Ablation: Pure OpenMP over Morton-sorted array
#pragma omp parallel for schedule(dynamic, 64) num_threads(config->n_threads)
    for (int i = 0; i < N; i++) {
      compute_force_single(sys->pos_x[i], sys->pos_y[i], sys->mass[i], i, root,
                           &sys->fx[i], &sys->fy[i], G, config->theta_max);
    }
  }

  if (config->use_arena) {
    // Arena reset happens at start of next call or above
  } else {
    free_tree_rec(root);
  }
}
