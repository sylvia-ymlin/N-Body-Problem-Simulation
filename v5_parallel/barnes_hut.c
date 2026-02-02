#include "barnes_hut.h"

/** Barnes-Hut algorithm
 * Based on the clustered particles, utilize the spatial and temporal locality
 * to improve the cache performance.
 * -------------------------------------------------- */
void barnes_hut(double *pos_x, double *pos_y, double *mass, int N, int *cluster,
                double *region, int *clusters_size, int k, double *fx,
                double *fy, int n_threads, double THETA_MAX, NodeArena *arena) {
  /* initialize the forces */
  for (int i = 0; i < N; i++) {
    fx[i] = 0;
    fy[i] = 0;
  }

  /* Reset the memory arena for the new tree build */
  reset_arena(arena);

  /* Build the global tree */
  /* Initialize the root with 2D region */
  TNode *tTree =
      create_new_TNode(arena, -1, region[0], region[1], region[2], region[3]);
  tTree->PID = 0;
  tTree->pos_x = pos_x[0];
  tTree->pos_y = pos_y[0];
  tTree->mass = mass[0];

  /* Insert each particle one by one */
  for (int i = 1; i < N; i++) {
    insert(arena, tTree, pos_x[i], pos_y[i], mass[i], i);
  }

  double G = 100.0 / N;
  /* Compute the forces */
  for (int i = 0; i < k; i++) {
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
    for (int j = 0; j < clusters_size[i]; j++) {
      int PID = cluster[i * N + j];
      compute_force_stackless(pos_x[PID], pos_y[PID], mass[PID], PID, tTree,
                              &fx[PID], &fy[PID], G, THETA_MAX);
    }
  }
}

/** Create a new tree node
 * According the position of the subsquare, determine the region of it.
 * -------------------------------------------------- */
TNode *create_new_TNode(NodeArena *arena, int index, double LB, double RB,
                        double DB, double UB) {
  if (arena->used >= arena->capacity) {
    fprintf(stderr, "Arena Overflow! Need capacity > %zu\n", arena->capacity);
    exit(1);
  }
  TNode *new_TNode = &arena->nodes[arena->used++];

  // Clear children pointers
  for (int i = 0; i < 4; i++)
    new_TNode->child[i] = NULL;
  new_TNode->PID = -1;
  new_TNode->mass = 0;

  double mid_x = 0.5 * (LB + RB);
  double mid_y = 0.5 * (DB + UB);

  new_TNode->LB = (index == 2 || index == 3) ? mid_x : LB;
  new_TNode->RB = (index == 2 || index == 3) ? RB : mid_x;
  new_TNode->DB = (index == 1 || index == 3) ? mid_y : DB;
  new_TNode->UB = (index == 1 || index == 3) ? UB : mid_y;

  if (index == -1) { // Root special case
    new_TNode->LB = LB; new_TNode->RB = RB;
    new_TNode->DB = DB; new_TNode->UB = UB;
  }
  return new_TNode;
}

/** Insert the particle into the tree */
int insert(NodeArena *arena, TNode *tNode, double pos_x, double pos_y,
           double mass, int PID) {
  double mid_x = 0.5 * (tNode->LB + tNode->RB);
  double mid_y = 0.5 * (tNode->DB + tNode->UB);
  double width = tNode->RB - tNode->LB;

  if (tNode->PID != -1) {
    if (width < 1e-12 || ((pos_x == tNode->pos_x) && (pos_y == tNode->pos_y))) {
        // Too small or identical position, just update mass and return
        double new_mass = tNode->mass + mass;
        tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
        tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
        tNode->mass = new_mass;
        return 0;
    }
    // Move existing particle to a child
    int old_index = (tNode->pos_y > mid_y) + 2 * (tNode->pos_x > mid_x);
    tNode->child[old_index] = create_new_TNode(
        arena, old_index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
    tNode->child[old_index]->PID = tNode->PID;
    tNode->child[old_index]->pos_x = tNode->pos_x;
    tNode->child[old_index]->pos_y = tNode->pos_y;
    tNode->child[old_index]->mass = tNode->mass;
    tNode->PID = -1;
  }

  /* Update CoM (Center of Mass) */
  double new_mass = tNode->mass + mass;
  tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
  tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
  tNode->mass = new_mass;

  int index = (pos_y > mid_y) + 2 * (pos_x > mid_x);
  if (tNode->child[index] == NULL) {
    tNode->child[index] =
        create_new_TNode(arena, index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
    tNode->child[index]->pos_x = pos_x;
    tNode->child[index]->pos_y = pos_y;
    tNode->child[index]->mass = mass;
    tNode->child[index]->PID = PID;
  } else {
    insert(arena, tNode->child[index], pos_x, pos_y, mass, PID);
  }
  return 0;
}

static inline void _compute_force(double pos_x, double pos_y, double mass,
                                  TNode *tNode, double *fx, double *fy,
                                  double G) {
  double r_x = pos_x - tNode->pos_x;
  double r_y = pos_y - tNode->pos_y;
  double r_sq = r_x * r_x + r_y * r_y + (EPSILON_O * EPSILON_O);
  double r_inv = 1.0 / sqrt(r_sq);
  double r_inv3 = r_inv * r_inv * r_inv;
  double force_factor = -G * mass * tNode->mass * r_inv3;
  *fx += force_factor * r_x;
  *fy += force_factor * r_y;
}


/** Compute the force between the particle and the node (Stackless) */
void compute_force_stackless(double pos_x, double pos_y, double mass, int PID,
                             TNode *root, double *fx, double *fy, double G,
                             double THETA_MAX) {
  TNode *stack[256];
  int sp = 0;

  if (root)
    stack[sp++] = root;

  while (sp > 0) {
    TNode *tNode = stack[--sp];

    if (tNode->PID == PID)
      continue;

    if (tNode->PID != -1) {
      _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
    } else {
      double mid_x = 0.5 * (tNode->LB + tNode->RB);
      double mid_y = 0.5 * (tNode->DB + tNode->UB);
      double width = (tNode->RB - tNode->LB);
      double dist_sq =
          (pos_x - mid_x) * (pos_x - mid_x) + (pos_y - mid_y) * (pos_y - mid_y);

      if (width * width <= THETA_MAX * THETA_MAX * dist_sq) {
        _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
      } else {
        for (int i = 3; i >= 0; i--) {
          if (tNode->child[i])
            stack[sp++] = tNode->child[i];
        }
      }
    }
  }
}

/* Arena Management Implementation */
void init_arena(NodeArena *arena, size_t capacity) {
  arena->nodes = (TNode *)calloc(capacity, sizeof(TNode));
  if (!arena->nodes) {
    fprintf(stderr, "Fatal: Failed to allocate arena with capacity %zu\n",
            capacity);
    exit(1);
  }
  arena->capacity = capacity;
  arena->used = 0;
}

void free_arena(NodeArena *arena) {
  if (arena->nodes) {
    free(arena->nodes);
    arena->nodes = NULL;
  }
  arena->capacity = 0;
  arena->used = 0;
}

void reset_arena(NodeArena *arena) {
  arena->used = 0;
}
