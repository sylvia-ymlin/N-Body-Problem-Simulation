#include "barnes_hut.h"
#include <math.h>
#include <string.h>

static TNode *create_node(double x_min, double x_max, double y_min,
                          double y_max) {
  TNode *node = (TNode *)malloc(sizeof(TNode));
  node->mass = 0;
  node->x = 0;
  node->y = 0;
  node->x_min = x_min;
  node->x_max = x_max;
  node->y_min = y_min;
  node->y_max = y_max;
  node->is_leaf = 1;
  node->particle_idx = -1;
  for (int i = 0; i < 4; i++)
    node->children[i] = NULL;
  return node;
}

static void free_tree(TNode *node) {
  if (!node)
    return;
  for (int i = 0; i < 4; i++)
    free_tree(node->children[i]);
  free(node);
}

static void insert(TNode *node, int idx, double *px, double *py, double *mass) {
  if (node->mass == 0 && node->is_leaf && node->particle_idx == -1) {
    node->mass = mass[idx];
    node->x = px[idx];
    node->y = py[idx];
    node->particle_idx = idx;
    return;
  }

  if (node->is_leaf) {
    int old_idx = node->particle_idx;
    node->is_leaf = 0;
    node->particle_idx = -1;
    if (old_idx != -1) {
      insert(node, old_idx, px, py, mass);
    }
  }

  double mid_x = (node->x_min + node->x_max) * 0.5;
  double mid_y = (node->y_min + node->y_max) * 0.5;
  int quad = 0;
  if (px[idx] > mid_x)
    quad += 1;
  if (py[idx] > mid_y)
    quad += 2;

  if (!node->children[quad]) {
    double nx_min = (quad & 1) ? mid_x : node->x_min;
    double nx_max = (quad & 1) ? node->x_max : mid_x;
    double ny_min = (quad & 2) ? mid_y : node->y_min;
    double ny_max = (quad & 2) ? node->y_max : mid_y;
    node->children[quad] = create_node(nx_min, nx_max, ny_min, ny_max);
  }
  insert(node->children[quad], idx, px, py, mass);

  // Update center of mass
  double total_mass = node->mass + mass[idx];
  node->x = (node->x * node->mass + px[idx] * mass[idx]) / total_mass;
  node->y = (node->y * node->mass + py[idx] * mass[idx]) / total_mass;
  node->mass = total_mass;
}

static void compute_force(TNode *node, int idx, double *px, double *py,
                          double *mass, double *fx, double *fy,
                          double THETA_MAX, double G) {
  if (!node || node->mass == 0 || node->particle_idx == idx)
    return;

  double dx = node->x - px[idx];
  double dy = node->y - py[idx];
  double r = sqrt(dx * dx + dy * dy + 1e-6);
  double s = node->x_max - node->x_min;

  if (node->is_leaf || (s / r < THETA_MAX)) {
    double f = G * mass[idx] * node->mass / (r * r * r);
    *fx += f * dx;
    *fy += f * dy;
  } else {
    for (int i = 0; i < 4; i++) {
      compute_force(node->children[i], idx, px, py, mass, fx, fy, THETA_MAX, G);
    }
  }
}

void barnes_hut(double *px, double *py, double *mass, int N, double *fx,
                double *fy, double THETA_MAX) {
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

  TNode *root = create_node(x_min, x_max, y_min, y_max);
  for (int i = 0; i < N; i++) {
    insert(root, i, px, py, mass);
  }

  double G = 100.0 / N;
  for (int i = 0; i < N; i++) {
    fx[i] = 0;
    fy[i] = 0;
    compute_force(root, i, px, py, mass, &fx[i], &fy[i], THETA_MAX, G);
  }

  free_tree(root);
}
