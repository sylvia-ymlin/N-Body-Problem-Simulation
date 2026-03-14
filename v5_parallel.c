#include "simulation.h"
#include "src/ds.h"
#include "src/kmeans.h"
#include "src/morton.h"
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>

/**
 * v5: Parallel Load Balancing
 *
 * Final architecture utilizing OpenMP dynamic scheduling over Morton-sorted
 * data. Spatial sequence provides cache-safe chunks for thread work-stealing.
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

static void compute_force_single(double pos_x, double pos_y, double mass,
                                 int PID, TNode *root, double *res_fx,
                                 double *res_fy, double G, double THETA) {
  TNode *stack[256];
  int sp = 0;
  if (root)
    stack[sp++] = root;

  double local_fx = 0;
  double local_fy = 0;

  while (sp > 0) {
    TNode *node = stack[--sp];
    if (node->PID == PID)
      continue;

    double dx = pos_x - node->pos_x;
    double dy = pos_y - node->pos_y;
    double r_sq = dx * dx + dy * dy;
    double s = node->RB - node->LB;

    if (node->PID != -1 || (s * s < THETA * THETA * r_sq)) {
      double r_inv = 1.0 / sqrt(r_sq + EPSILON * EPSILON);
      double f = G * mass * node->mass * (r_inv * r_inv * r_inv);
      local_fx += f * (-dx);
      local_fy += f * (-dy);
    } else {
      for (int i = 0; i < 4; i++) {
        if (node->child[i])
          stack[sp++] = node->child[i];
      }
    }
  }

  // Final single write reduces cache-line contention
  *res_fx = local_fx;
  *res_fy = local_fy;
}

void compute_force_v5_parallel(ParticleSystem *sys, KernelConfig *config) {
  int N = sys->N;
  int k = config->k_clusters;

  // 1. Build Bounding Box
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
  double dx = x_max - x_min;
  double dy = y_max - y_min;
  double d = (dx > dy) ? dx : dy;
  x_max = x_min + d;
  y_max = y_min + d;
  x_min -= d * 0.05;
  x_max += d * 0.05;
  y_min -= d * 0.05;
  y_max += d * 0.05;

  // 2. Parallel Strategy Selection (Ablation Study)
  if (config->k_clusters <= 0) {
    // Optimal for shared memory: Morton Sort + Dynamic Scheduling
    z_order_sort(sys, x_min, x_max, y_min, y_max);
  } else {
    // K-Means Geometric Partitioning
    static int *clusters = NULL;
    static int *c_size = NULL;
    static int last_N = 0, last_k = 0;
    if (!clusters || last_N < N) {
      free(clusters);
      clusters = (int *)malloc(N * sizeof(int));
      last_N = N;
    }
    if (!c_size || last_k < config->k_clusters) {
      free(c_size);
      c_size = (int *)malloc(config->k_clusters * sizeof(int));
      last_k = config->k_clusters;
    }
    kmeans(sys, clusters, c_size, config->k_clusters, config->n_threads);
  }

  // 3. Build Tree (Serial - Bottleneck but required for correctness)
  double t_tree_start = get_sim_time();
  if (arena.buffer == NULL || arena_N != N) {
    if (arena.buffer)
      free_arena(&arena);
    init_arena(&arena, N * 100);
    arena_N = N;
  }
  reset_arena(&arena);
  TNode *root = create_node(x_min, x_max, y_min, y_max, &arena);
  for (int i = 0; i < N; i++)
    insert(root, i, sys, &arena);
  double t_tree_end = get_sim_time();

  // 4. Compute Forces (Parallel)
  double t_force_start = get_sim_time();
  double G = 100.0 / N;
#pragma omp parallel for schedule(dynamic, 128) num_threads(config->n_threads)
  for (int i = 0; i < N; i++) {
    compute_force_single(sys->pos_x[i], sys->pos_y[i], sys->mass[i], i, root,
                         &sys->fx[i], &sys->fy[i], G, config->theta_max);
  }
  double t_force_end = get_sim_time();

  // Report hotspots (visible in bench Mode)
  static int first = 1;
  if (first) {
    printf("\n[HP] Tree Build: %.4f s | Force Compute: %.4f s\n",
           t_tree_end - t_tree_start, t_force_end - t_force_start);
    // first = 0; // Show every step for now
  }
}
