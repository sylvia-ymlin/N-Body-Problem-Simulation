#include "ds.h"
#include "kmeans.h"
#include "morton.h"
#include "types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define G_FACTOR 100.0
#define EPSILON  1e-3

static const double COINCIDENT_EPS      = 1e-9;
static const double DOMAIN_PADDING_FRAC = 0.05;
static const int    ARENA_NODE_FACTOR   = 100;
#define CHUNK_SIZE 128

/* Pre-allocated node pool, reused every timestep */
static NodeArena arena   = {NULL, 0, 0};
static int       arena_N = 0;

/* Shared within each timestep */
static double G_val     = 0.0;
static double theta_val = 0.0;

static int    is_leaf(TNode* node);
static int    quadrant(double px, double py, double mx, double my);
static void   init_domain(const double* x, const double* y, int N,
                          double* x_min, double* x_max,
                          double* y_min, double* y_max);
static void   expand_domain_if_needed(const double* x, const double* y, int N,
                                      double* x_min, double* x_max,
                                      double* y_min, double* y_max);
static TNode* create_node(double LB, double RB, double DB, double UB);
static void   insert(TNode* node, int idx, ParticleSystem* sys);
static void   compute_force_single(int i, ParticleSystem* sys, TNode* root,
                                   double* res_fx, double* res_fy);

void compute_force_barnes_hut(ParticleSystem* sys, KernelConfig* config) {
    int N = sys->N;
    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    double* fx_out  = sys->fx;
    double* fy_out  = sys->fy;

    static int    domain_initialized = 0;
    static int    domain_N           = 0;
    static double x_min = 0.0, x_max = 0.0, y_min = 0.0, y_max = 0.0;

    if (!domain_initialized || domain_N != N) {
        init_domain(x, y, N, &x_min, &x_max, &y_min, &y_max);
        domain_initialized = 1;
        domain_N = N;
    } else {
        expand_domain_if_needed(x, y, N, &x_min, &x_max, &y_min, &y_max);
    }

    /* Reorder particles for better cache locality during tree traversal */
    if (config->k_clusters > 0) {
        static int*   clusters            = NULL;
        static int*   c_size              = NULL;
        static int    last_N              = 0, last_k = 0;
        static double last_recluster_time = -1.0;
        static const double RECLUSTER_INTERVAL = 1e-4;

        if (!clusters || last_N < N) {
            free(clusters);
            clusters = (int*)malloc(N * sizeof(int));
            last_N   = N;
        }
        if (!c_size || last_k < config->k_clusters) {
            free(c_size);
            c_size = (int*)malloc(config->k_clusters * sizeof(int));
            last_k = config->k_clusters;
        }
        if (last_recluster_time < 0.0 ||
            config->current_time - last_recluster_time >= RECLUSTER_INTERVAL) {
            kmeans(sys, clusters, c_size, config->k_clusters, config->n_threads);
            last_recluster_time = config->current_time;
        }
    } else {
        z_order_sort(sys, x_min, x_max, y_min, y_max);
    }

    /* Allocate (or resize) the arena once; reset it each timestep */
    if (arena.buffer == NULL || arena_N != N) {
        if (arena.buffer) free_arena(&arena);
        init_arena(&arena, N * ARENA_NODE_FACTOR);
        arena_N = N;
    }
    reset_arena(&arena);

    G_val     = G_FACTOR / N;
    theta_val = config->theta_max;

    TNode* root = create_node(x_min, x_max, y_min, y_max);
    for (int i = 0; i < N; i++)
        insert(root, i, sys);

    /* Only the per-particle force traversal is parallelised;
     * tree building and Morton sort remain serial. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE) num_threads(config->n_threads)
#endif
    for (int i = 0; i < N; i++)
        compute_force_single(i, sys, root, &fx_out[i], &fy_out[i]);
}

/** Compute the bounding square for all particles with a small padding.
 * ----------------------------------------------------------------- */
static void init_domain(const double* x, const double* y, int N,
                        double* x_min, double* x_max,
                        double* y_min, double* y_max) {
    double nx_min = x[0], nx_max = x[0];
    double ny_min = y[0], ny_max = y[0];
    for (int i = 1; i < N; i++) {
        if (x[i] < nx_min) nx_min = x[i];
        if (x[i] > nx_max) nx_max = x[i];
        if (y[i] < ny_min) ny_min = y[i];
        if (y[i] > ny_max) ny_max = y[i];
    }

    double size = nx_max - nx_min;
    if ((ny_max - ny_min) > size) size = ny_max - ny_min;
    if (size <= 0.0) size = 1e-6;

    *x_min = nx_min - size * DOMAIN_PADDING_FRAC;
    *x_max = nx_min + size + size * DOMAIN_PADDING_FRAC;
    *y_min = ny_min - size * DOMAIN_PADDING_FRAC;
    *y_max = ny_min + size + size * DOMAIN_PADDING_FRAC;
}

/** Grow the domain if any particle has escaped the current bounds.
 * ----------------------------------------------------------------- */
static void expand_domain_if_needed(const double* x, const double* y, int N,
                                    double* x_min, double* x_max,
                                    double* y_min, double* y_max) {
    for (int i = 0; i < N; i++) {
        if (x[i] < *x_min || x[i] > *x_max || y[i] < *y_min || y[i] > *y_max) {
            init_domain(x, y, N, x_min, x_max, y_min, y_max);
            return;
        }
    }
}

/** Iterative tree traversal using an explicit stack.
 * Accepts a subtree as one pseudo-body when cell_width / r < theta.
 * ----------------------------------------------------------------- */
static void compute_force_single(int i, ParticleSystem* sys, TNode* root,
                                 double* res_fx, double* res_fy) {
    double pos_x = sys->pos_x[i];
    double pos_y = sys->pos_y[i];
    double mass  = sys->mass[i];

    TNode* stack[256];
    int sp = 0;
    if (root) stack[sp++] = root;

    double fx = 0.0, fy = 0.0;

    while (sp > 0) {
        TNode* node = stack[--sp];
        if (node->particle_idx == i) continue;

        double dx = pos_x - node->pos_x;
        double dy = pos_y - node->pos_y;
        double r  = sqrt(dx * dx + dy * dy);
        double s  = node->x_max - node->x_min;

        if (node->particle_idx != -1 || (s < theta_val * r)) {
            double denom = r + EPSILON;
            double f = G_val * mass * node->mass / (denom * denom * denom);
            fx += f * (-dx);
            fy += f * (-dy);
        } else {
            for (int j = 0; j < 4; j++) {
                if (node->child[j]) stack[sp++] = node->child[j];
            }
        }
    }

    *res_fx = fx;
    *res_fy = fy;
}

/** Insert particle idx into the quadtree.
 * Coincident particles are merged into one aggregate leaf.
 * ----------------------------------------------------------------- */
static void insert(TNode* node, int idx, ParticleSystem* sys) {
    double px   = sys->pos_x[idx];
    double py   = sys->pos_y[idx];
    double mass = sys->mass[idx];

    /* Empty leaf: just store this particle */
    if (node->particle_idx == -1 && is_leaf(node) && node->mass == 0) {
        node->particle_idx = idx;
        node->mass  = mass;
        node->pos_x = px;
        node->pos_y = py;
        return;
    }

    if (is_leaf(node)) {
        int old = node->particle_idx;
        if (old != -1) {
            /* Merge particles that land on the same spot to avoid infinite splits */
            if (fabs(px - sys->pos_x[old]) < COINCIDENT_EPS &&
                fabs(py - sys->pos_y[old]) < COINCIDENT_EPS) {
                double mt = node->mass + mass;
                node->pos_x = (node->pos_x * node->mass + px * mass) / mt;
                node->pos_y = (node->pos_y * node->mass + py * mass) / mt;
                node->mass  = mt;
                return;
            }

            /* Subdivide: push the existing particle down into a child */
            node->particle_idx = -1;
            double old_px = node->pos_x, old_py = node->pos_y, old_m = node->mass;
            node->mass  = 0;
            node->pos_x = 0;
            node->pos_y = 0;

            double mx = (node->x_min + node->x_max) * 0.5;
            double my = (node->y_min + node->y_max) * 0.5;
            int oq = quadrant(old_px, old_py, mx, my);
            double lb, rb, db, ub;
            if (oq == 0)      { lb = node->x_min; rb = mx;          db = node->y_min; ub = my;          }
            else if (oq == 1) { lb = mx;          rb = node->x_max; db = node->y_min; ub = my;          }
            else if (oq == 2) { lb = node->x_min; rb = mx;          db = my;          ub = node->y_max; }
            else              { lb = mx;          rb = node->x_max; db = my;          ub = node->y_max; }
            node->child[oq] = create_node(lb, rb, db, ub);
            insert(node->child[oq], old, sys);

            node->mass  = old_m;
            node->pos_x = old_px;
            node->pos_y = old_py;
        }
    }

    /* Insert new particle into the appropriate child cell */
    double mx = (node->x_min + node->x_max) * 0.5;
    double my = (node->y_min + node->y_max) * 0.5;
    int q = quadrant(px, py, mx, my);
    if (!node->child[q]) {
        double lb, rb, db, ub;
        if (q == 0)      { lb = node->x_min; rb = mx;          db = node->y_min; ub = my;          }
        else if (q == 1) { lb = mx;          rb = node->x_max; db = node->y_min; ub = my;          }
        else if (q == 2) { lb = node->x_min; rb = mx;          db = my;          ub = node->y_max; }
        else             { lb = mx;          rb = node->x_max; db = my;          ub = node->y_max; }
        node->child[q] = create_node(lb, rb, db, ub);
    }
    insert(node->child[q], idx, sys);

    /* Update this internal node's mass and centre of mass */
    double mt = node->mass + mass;
    node->pos_x = (node->pos_x * node->mass + px * mass) / mt;
    node->pos_y = (node->pos_y * node->mass + py * mass) / mt;
    node->mass  = mt;
}

/* Allocate a node from the arena instead of calling malloc each time */
static TNode* create_node(double LB, double RB, double DB, double UB) {
    TNode* node = arena_alloc(&arena);
    if (!node) {
        fprintf(stderr, "Error: arena out of memory!\n");
        exit(1);
    }
    node->x_min = LB;  node->x_max = RB;
    node->y_min = DB;  node->y_max = UB;
    node->pos_x = 0;   node->pos_y = 0;
    node->mass  = 0;
    node->particle_idx = -1;
    for (int i = 0; i < 4; i++)
        node->child[i] = NULL;
    return node;
}

static int is_leaf(TNode* node) {
    return node->child[0] == NULL && node->child[1] == NULL &&
           node->child[2] == NULL && node->child[3] == NULL;
}

/** Determine which quadrant (0=SW 1=SE 2=NW 3=NE) a point falls in.
 * ----------------------------------------------------------------- */
static int quadrant(double px, double py, double mx, double my) {
    if (px <= mx && py <= my) return 0; /* SW */
    if (px >  mx && py <= my) return 1; /* SE */
    if (px <= mx && py >  my) return 2; /* NW */
    return 3;                            /* NE */
}
