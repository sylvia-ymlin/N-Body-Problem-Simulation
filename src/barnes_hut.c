#include "barnes_hut.h"

/** Barnes-Hut algorithm
 * Based on the clustered particles, utilize the spatial and temporal locality
 * to improve the cache performance.
 * - Build the local tree for each cluster, parallelize the insertion of the
 * particles => no significant improvement.
 * - Or build the global tree
 *  -------------------------------------------------- */
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
  /* Initialize the root */
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

  /* Compute the forces */
  for (int i = 0; i < k; i++) {
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(n_threads)
    for (int j = 0; j < clusters_size[i]; j++) {
      /* each particle traverses all the trees */
      int PID = cluster[i * N + j];
      compute_force_stackless(pos_x[PID], pos_y[PID], mass[PID], PID, tTree,
                              &fx[PID], &fy[PID], N, THETA_MAX);
    }
  }

  /* Tree is destroyed automatically when arena is reused/reset next time */
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

  // Clear children pointers explicitly since memory is reused
  new_TNode->child[0] = NULL;
  new_TNode->child[1] = NULL;
  new_TNode->child[2] = NULL;
  new_TNode->child[3] = NULL;
  new_TNode->PID = -1;
  new_TNode->mass = 0;

  double mid_x = 0.5 * (LB + RB);
  double mid_y = 0.5 * (DB + UB);
  switch (index) {
  case -1:
    new_TNode->LB = LB;
    new_TNode->RB = RB;
    new_TNode->DB = DB;
    new_TNode->UB = UB;
    break;
  case 0:
    new_TNode->LB = LB;
    new_TNode->RB = mid_x;
    new_TNode->DB = DB;
    new_TNode->UB = mid_y;
    break;
  case 1:
    new_TNode->LB = LB;
    new_TNode->RB = mid_x;
    new_TNode->DB = mid_y;
    new_TNode->UB = UB;
    break;
  case 2:
    new_TNode->LB = mid_x;
    new_TNode->RB = RB;
    new_TNode->DB = DB;
    new_TNode->UB = mid_y;
    break;
  case 3:
    new_TNode->LB = mid_x;
    new_TNode->RB = RB;
    new_TNode->DB = mid_y;
    new_TNode->UB = UB;
    break;
  }
  return new_TNode;
}

/** Insert the particle into the tree
 * If the node represents a particle, then split the node and insert the
 * particle into the corresponding child node. If the node is representing a big
 * object, then update the massive and the centroid of the node.
 * - determine the region the particle belongs to and insert it into the
 * corresponding child node.
 * - if the child node is NULL, then create a new node and insert the particle
 * into it.
 * - else, recursively insert the particle into the child node.
 * -------------------------------------------------- */
int insert(NodeArena *arena, TNode *tNode, double pos_x, double pos_y,
           double mass, int PID) {
  double mid_x = 0.5 * (tNode->LB + tNode->RB);
  double mid_y = 0.5 * (tNode->UB + tNode->DB);

  if (tNode->PID != -1) {
    if ((pos_x == tNode->pos_x) && (pos_y == tNode->pos_y)) {
      // Critical collision detected. For stability, we terminate.
      printf("Two particles are detected at the same location and the "
             "simulation terminates.\n");
      exit(0);
    }
    int index = (tNode->pos_y > mid_y) + 2 * (tNode->pos_x > mid_x);
    tNode->child[index] = create_new_TNode(arena, index, tNode->LB, tNode->RB,
                                           tNode->DB, tNode->UB);
    tNode->child[index]->PID = tNode->PID;
    tNode->child[index]->pos_x = tNode->pos_x;
    tNode->child[index]->pos_y = tNode->pos_y;
    tNode->child[index]->mass = tNode->mass;
    tNode->PID = -1;
  }

  /* Update the massive and the centroid */
  double new_mass = tNode->mass + mass;
  tNode->pos_x = (mass * pos_x + tNode->mass * tNode->pos_x) / new_mass;
  tNode->pos_y = (mass * pos_y + tNode->mass * tNode->pos_y) / new_mass;
  tNode->mass = new_mass;

  int index = (pos_y > mid_y) + 2 * (pos_x > mid_x);
  if (tNode->child[index] == NULL) {
    tNode->child[index] = create_new_TNode(arena, index, tNode->LB, tNode->RB,
                                           tNode->DB, tNode->UB);
    tNode->child[index]->pos_x = pos_x;
    tNode->child[index]->pos_y = pos_y;
    tNode->child[index]->mass = mass;
    tNode->child[index]->PID = PID;
  } else {
    insert(arena, tNode->child[index], pos_x, pos_y, mass, PID);
  }
  return 0;
}

void _compute_force(double pos_x, double pos_y, double mass, TNode *tNode,
                    double *fx, double *fy, double G) {
  double r_x = pos_x - tNode->pos_x;
  double r_y = pos_y - tNode->pos_y;
  double r_squared = r_x * r_x + r_y * r_y;
  double r_plummer = sqrt(r_squared) + EPSILON_O;
  double force_factor =
      -G * mass * tNode->mass / (r_plummer * r_plummer * r_plummer);
  *fx += force_factor * r_x;
  *fy += force_factor * r_y;
}

/** Compute the force between the particle and the node
 * If the node is a leaf node, directly calculate the force.
 * If the node is a big object, calculate the theta to check if the node is far
 * enough, then calculate the force.
 * - calculate the distance between the particle and the node.
 * - calculate the force factor and update the force.
 * -------------------------------------------------- */
/** Compute the force between the particle and the node (Recursive) */
void compute_force_recursive(double pos_x, double pos_y, double mass, int PID,
                             TNode *tNode, double *fx, double *fy, int N,
                             double THETA_MAX) {
  double G = 100.0 / N;
  if (!tNode || tNode->PID == PID)
    return;

  if (tNode->PID != -1) {
    _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
  } else {
    double mid_x = 0.5 * (tNode->LB + tNode->RB);
    double mid_y = 0.5 * (tNode->DB + tNode->UB);
    double width = (tNode->RB - tNode->LB);
    double theta = width / sqrt((pos_x - mid_x) * (pos_x - mid_x) +
                                (pos_y - mid_y) * (pos_y - mid_y));

    if (theta <= THETA_MAX) {
      _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
    } else {
      compute_force_recursive(pos_x, pos_y, mass, PID, tNode->child[0], fx, fy,
                              N, THETA_MAX);
      compute_force_recursive(pos_x, pos_y, mass, PID, tNode->child[1], fx, fy,
                              N, THETA_MAX);
      compute_force_recursive(pos_x, pos_y, mass, PID, tNode->child[2], fx, fy,
                              N, THETA_MAX);
      compute_force_recursive(pos_x, pos_y, mass, PID, tNode->child[3], fx, fy,
                              N, THETA_MAX);
    }
  }
}

/** Compute the force between the particle and the node (Stackless) */
void compute_force_stackless(double pos_x, double pos_y, double mass, int PID,
                             TNode *root, double *fx, double *fy, int N,
                             double THETA_MAX) {
  double G = 100.0 / N;

  // Explicit Stack
  TNode *stack[256];
  int sp = 0;

  // Push root
  if (root)
    stack[sp++] = root;

  while (sp > 0) {
    TNode *tNode = stack[--sp]; // Pop

    if (tNode->PID == PID)
      continue;

    if (tNode->PID != -1) { // Leaf node (Particle)
      _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
    } else { // Internal node
      double mid_x = 0.5 * (tNode->LB + tNode->RB);
      double mid_y = 0.5 * (tNode->DB + tNode->UB);
      double width = (tNode->RB - tNode->LB);
      double dist_sq =
          (pos_x - mid_x) * (pos_x - mid_x) + (pos_y - mid_y) * (pos_y - mid_y);
      double theta = width / sqrt(dist_sq);

      if (theta <= THETA_MAX) {
        _compute_force(pos_x, pos_y, mass, tNode, fx, fy, G);
      } else {
        // Push children (order doesn't strictly matter for correctness,
        // but reverse order simulates recursion order: 3,2,1,0 -> pop 0 first)
        if (tNode->child[3])
          stack[sp++] = tNode->child[3];
        if (tNode->child[2])
          stack[sp++] = tNode->child[2];
        if (tNode->child[1])
          stack[sp++] = tNode->child[1];
        if (tNode->child[0])
          stack[sp++] = tNode->child[0];
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
  // Nodes are overwritten as needed, no need to zero out everything for
  // performance
}