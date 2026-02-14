#include "barnes_hut.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CHUNK_SIZE 8

static void free_tree(TNode *node) {
  if (!node) return;
  for (int i = 0; i < 4; i++) {
    if (node->child[i]) {
      free_tree(node->child[i]);
    }
  }
  free(node);
}

// Create a new tree node representing a sub-quadrant
TNode *create_new_TNode(int index, double LB, double RB,
                        double DB, double UB, NodeArena *arena) {
  TNode *new_TNode;
  if (arena) {
    new_TNode = arena_alloc(arena);
  } else {
    new_TNode = (TNode *)malloc(sizeof(TNode));
  }

  if (!new_TNode) {
      fprintf(stderr, "Memory allocation failed!\n");
      exit(1);
  }

  // Clear children pointers
  for (int i = 0; i < 4; i++)
    new_TNode->child[i] = NULL;
  new_TNode->PID = -1;
  new_TNode->mass = 0;
  new_TNode->pos_x = 0;
  new_TNode->pos_y = 0;

  double mid_x = 0.5 * (LB + RB);
  double mid_y = 0.5 * (DB + UB);

  if (index == -1) { // Root special case
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

/** Insert the particle into the tree */
int insert(TNode *tNode, double pos_x, double pos_y,
           double mass, int PID, NodeArena *arena) {
  double mid_x = 0.5 * (tNode->LB + tNode->RB);
  double mid_y = 0.5 * (tNode->DB + tNode->UB);
  double width = tNode->RB - tNode->LB;

  // If we are at a leaf node that is occupied
  if (tNode->PID != -1) {
    // If the node is very small or particles are at the exact same position
    if (width < 1e-12 || ((pos_x == tNode->pos_x) && (pos_y == tNode->pos_y))) {
      // Just update mass and return (collision/overlap handling)
      double new_mass = tNode->mass + mass;
      tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
      tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
      tNode->mass = new_mass;
      return 0;
    }
    
    // Otherwise, push the existing particle down to a child
    int old_index = (tNode->pos_y > mid_y) + 2 * (tNode->pos_x > mid_x);
    tNode->child[old_index] = create_new_TNode(old_index, tNode->LB,
                                               tNode->RB, tNode->DB, tNode->UB, arena);
    tNode->child[old_index]->PID = tNode->PID;
    tNode->child[old_index]->pos_x = tNode->pos_x;
    tNode->child[old_index]->pos_y = tNode->pos_y;
    tNode->child[old_index]->mass = tNode->mass;
    
    // Mark current node as internal
    tNode->PID = -1;
  }

  /* Update CoM (Center of Mass) for this node */
  double new_mass = tNode->mass + mass;
  // If this is a fresh node (mass=0), avoid div by zero if mass is also 0 (unlikely)
  if (new_mass > 0) {
      tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
      tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
  } else {
      tNode->pos_x = pos_x;
      tNode->pos_y = pos_y;
  }
  tNode->mass = new_mass;

  // Determine which child the new particle belongs to
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
  return 0;
}

static inline void _compute_force(double pos_x, double pos_y, double mass,
                                  TNode *tNode, double *fx, double *fy,
                                  double G) {
  double r_x = tNode->pos_x - pos_x;
  double r_y = tNode->pos_y - pos_y;
  double r_sq = r_x * r_x + r_y * r_y + (EPSILON_O * EPSILON_O);
  double r_inv = 1.0 / sqrt(r_sq);
  double r_inv3 = r_inv * r_inv * r_inv;
  double force_factor = G * mass * tNode->mass * r_inv3; 
  *fx += force_factor * r_x;
  *fy += force_factor * r_y;
}

// Compute gravitational force using the BH tree
void compute_force(double pos_x, double pos_y, double mass, int PID,
                   TNode *root, double *fx, double *fy, double G,
                   double THETA_MAX) {
  // Non-recursive stack-based traversal to avoid recursion depth issues
  TNode *stack[2048]; // Increased stack size
  int sp = 0;

  if (root)
    stack[sp++] = root;

  while (sp > 0) {
    TNode *tNode = stack[--sp];

    // Skip if it's the same particle
    if (tNode->PID == PID)
      continue;

    double mid_x = 0.5 * (tNode->LB + tNode->RB);
    double mid_y = 0.5 * (tNode->DB + tNode->UB);
    double width = (tNode->RB - tNode->LB);
    double dx = pos_x - tNode->pos_x;
    double dy = pos_y - tNode->pos_y;
    double dist_sq = dx*dx + dy*dy;

    // If it's a leaf node (PID != -1), compute force directly
    if (tNode->PID != -1) {
      _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
    } else {
      // It's an internal node
      // Check Barnes-Hut condition: s/d < theta
      // s = width, d = sqrt(dist_sq)
      // s^2 / d^2 < theta^2  => s^2 < theta^2 * d^2
      
      if (width * width < THETA_MAX * THETA_MAX * dist_sq) {
        // Far enough, treat as point mass
        _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
      } else {
        // Too close, open the node (push children)
        for (int i = 0; i < 4; i++) {
          if (tNode->child[i]) {
            if (sp < 2048) {
              stack[sp++] = tNode->child[i];
            } else {
              // Fallback or error? For now, just ignore to avoid crash, but warn
              // fprintf(stderr, "Stack overflow\n");
            }
          }
        }
      }
    }
  }
}

void barnes_hut(double *pos_x, double *pos_y, double *mass, int N, int *cluster,
                double *region, int *clusters_size, int k, double *fx,
                double *fy, int n_threads, double THETA_MAX, int use_arena) {
  /* initialize the forces */
  for (int i = 0; i < N; i++) {
    fx[i] = 0;
    fy[i] = 0;
  }

  /* Update region bounds dynamically */
  region[0] = pos_x[0];
  region[1] = pos_x[0];
  region[2] = pos_y[0];
  region[3] = pos_y[0];
  for (int i = 1; i < N; i++) {
    if (pos_x[i] < region[0])
      region[0] = pos_x[i];
    if (pos_x[i] > region[1])
      region[1] = pos_x[i];
    if (pos_y[i] < region[2])
      region[2] = pos_y[i];
    if (pos_y[i] > region[3])
      region[3] = pos_y[i];
  }
  double r_dx = region[1] - region[0];
  double r_dy = region[3] - region[2];
  if (r_dx < 1e-6)
    r_dx = 1e-6;
  if (r_dy < 1e-6)
    r_dy = 1e-6;
  region[0] -= r_dx * 0.05;
  region[1] += r_dx * 0.05;
  region[2] -= r_dy * 0.05;
  region[3] += r_dy * 0.05;

  NodeArena arena_struct;
  NodeArena *arena = NULL;
  
  if (use_arena) {
    // Heuristic: 2N nodes is usually enough for Barnes-Hut
    // But let's be safe with 4N or even more to avoid realloc (we don't have realloc in simple arena)
    init_arena(&arena_struct, N * 10); // Safe upper bound
    arena = &arena_struct;
  }

  /* Build the global tree */
  TNode *tTree = create_new_TNode(-1, region[0], region[1], region[2], region[3], arena);
  tTree->PID = 0;
  tTree->pos_x = pos_x[0];
  tTree->pos_y = pos_y[0];
  tTree->mass = mass[0];

  /* Insert each particle one by one */
  for (int i = 1; i < N; i++) {
    insert(tTree, pos_x[i], pos_y[i], mass[i], i, arena);
  }

  double G = 100.0 / N;
  
  /* Compute the forces */
  // Reconstruct offsets for packed clusters array
  int *offsets = (int *)malloc(k * sizeof(int));
  offsets[0] = 0;
  for (int i = 1; i < k; i++)
    offsets[i] = offsets[i - 1] + clusters_size[i - 1];

  #pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < clusters_size[i]; j++) {
      int PID = cluster[offsets[i] + j];
      compute_force(pos_x[PID], pos_y[PID], mass[PID], PID, tTree, &fx[PID],
                    &fy[PID], G, THETA_MAX);
    }
  }
  free(offsets);
  
  // IMPORTANT: Free the tree to fix memory leak
  if (use_arena) {
    free_arena(arena);
  } else {
    free_tree(tTree);
  }
}
