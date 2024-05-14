#include "core/domain.h"
#include "core/ds.h"
#include "core/kmeans.h"
#include "core/morton.h"
#include "simulation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static NodeArena arena = {NULL, 0, 0};
static int arena_N = 0;
static const int ARENA_NODE_FACTOR = 100;
static const double DOMAIN_PADDING_FRAC = 0.05;
static const double COINCIDENT_EPS = 1e-9;

static inline int is_leaf(TNode* node);
static inline int quadrant(double px, double py, double mx, double my);
static TNode* create_node(double LB, double RB, double DB, double UB,
                          NodeArena* a);
static void insert(TNode* node, int idx, ParticleSystem* sys, NodeArena* a);
static void compute_force_single(double pos_x, double pos_y, double mass,
                                 int particle_idx, TNode* root, double* res_fx,
                                 double* res_fy, double G, double THETA);

void compute_force_v4_morton(ParticleSystem* sys, KernelConfig* config) {
    int N = sys->N;
    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* m = sys->mass;
    double* fx_out = sys->fx;
    double* fy_out = sys->fy;

    static Domain domain = {0, 0, 0.0, 0.0, 0.0, 0.0};
    if (!domain.initialized || domain.N != N) {
        domain_square_init_from_particles(x, y, N, DOMAIN_PADDING_FRAC,
                                          &domain.x_min, &domain.x_max,
                                          &domain.y_min, &domain.y_max);
        domain.initialized = 1;
        domain.N = N;
    } else {
        domain_square_expand_if_needed(x, y, N, DOMAIN_PADDING_FRAC,
                                       &domain.x_min, &domain.x_max,
                                       &domain.y_min, &domain.y_max);
    }

    if (config->k_clusters <= 0) {
        z_order_sort(sys, domain.x_min, domain.x_max, domain.y_min, domain.y_max);
    } else {
        static int* clusters = NULL;
        static int* c_size = NULL;
        static int last_N = 0, last_k = 0;
        static double last_recluster_time = -1.0;
        static const double RECLUSTER_INTERVAL = 1e-4;

        if (!clusters || last_N < N) {
            free(clusters);
            clusters = (int*)malloc(N * sizeof(int));
            last_N = N;
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
    }

    if (arena.buffer == NULL || arena_N != sys->N) {
        if (arena.buffer)
            free_arena(&arena);
        init_arena(&arena, sys->N * ARENA_NODE_FACTOR);
        arena_N = sys->N;
    }
    reset_arena(&arena);

    TNode* root = create_node(domain.x_min, domain.x_max, domain.y_min,
                              domain.y_max, &arena);
    for (int i = 0; i < N; i++)
        insert(root, i, sys, &arena);

    double G = G_FACTOR / N;
    double theta = config->theta_max;

    for (int i = 0; i < N; i++) {
        compute_force_single(x[i], y[i], m[i], i, root, &fx_out[i], &fy_out[i],
                             G, theta);
    }
}

static inline int is_leaf(TNode* node) {
    return node->child[0] == NULL && node->child[1] == NULL &&
           node->child[2] == NULL && node->child[3] == NULL;
}

static inline int quadrant(double px, double py, double mx, double my) {
    return (px > mx) + 2 * (py > my);
}

static TNode* create_node(double LB, double RB, double DB, double UB,
                          NodeArena* a) {
    TNode* node = arena_alloc(a);
    if (!node) {
        fprintf(stderr, "Error: Arena out of memory!\n");
        exit(1);
    }
    node->x_min = LB;
    node->x_max = RB;
    node->y_min = DB;
    node->y_max = UB;
    node->pos_x = 0;
    node->pos_y = 0;
    node->mass = 0;
    node->particle_idx = -1;
    for (int i = 0; i < 4; i++)
        node->child[i] = NULL;
    return node;
}

static void insert(TNode* node, int idx, ParticleSystem* sys, NodeArena* a) {
    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* m = sys->mass;
    double px = x[idx], py = y[idx], mass = m[idx];

    if (node->particle_idx == -1 && is_leaf(node) && node->mass == 0) {
        node->particle_idx = idx;
        node->mass = mass;
        node->pos_x = px;
        node->pos_y = py;
        return;
    }

    if (is_leaf(node)) {
        int old = node->particle_idx;
        if (old != -1) {
            if (fabs(px - x[old]) < COINCIDENT_EPS &&
                fabs(py - y[old]) < COINCIDENT_EPS) {
                double mt = node->mass + mass;
                node->pos_x = (node->pos_x * node->mass + px * mass) / mt;
                node->pos_y = (node->pos_y * node->mass + py * mass) / mt;
                node->mass = mt;
                return;
            }
            node->particle_idx = -1;
            double old_px = node->pos_x, old_py = node->pos_y,
                   old_m = node->mass;
            node->mass = 0;
            node->pos_x = 0;
            node->pos_y = 0;

            double mx = (node->x_min + node->x_max) * 0.5,
                   my = (node->y_min + node->y_max) * 0.5;
            int oq = quadrant(old_px, old_py, mx, my);
            node->child[oq] = create_node(
                (oq & 1) ? mx : node->x_min, (oq & 1) ? node->x_max : mx,
                (oq & 2) ? my : node->y_min, (oq & 2) ? node->y_max : my, a);
            insert(node->child[oq], old, sys, a);

            node->mass = old_m;
            node->pos_x = old_px;
            node->pos_y = old_py;
        }
    }

    double mx = (node->x_min + node->x_max) * 0.5,
           my = (node->y_min + node->y_max) * 0.5;
    int q = quadrant(px, py, mx, my);
    if (!node->child[q]) {
        node->child[q] = create_node(
            (q & 1) ? mx : node->x_min, (q & 1) ? node->x_max : mx,
            (q & 2) ? my : node->y_min, (q & 2) ? node->y_max : my, a);
    }
    insert(node->child[q], idx, sys, a);

    double mt = node->mass + mass;
    node->pos_x = (node->pos_x * node->mass + px * mass) / mt;
    node->pos_y = (node->pos_y * node->mass + py * mass) / mt;
    node->mass = mt;
}

static void compute_force_single(double pos_x, double pos_y, double mass,
                                 int particle_idx, TNode* root, double* res_fx,
                                 double* res_fy, double G, double THETA) {
    TNode* stack[256];
    int sp = 0;
    if (root)
        stack[sp++] = root;

    double local_fx = 0;
    double local_fy = 0;

    while (sp > 0) {
        TNode* node = stack[--sp];
        if (node->particle_idx == particle_idx)
            continue;

        double dx = pos_x - node->pos_x;
        double dy = pos_y - node->pos_y;
        double r = sqrt(dx * dx + dy * dy);
        double s = node->x_max - node->x_min;

        const int is_single_particle = (node->particle_idx != -1);
        const int accept_approx =
            is_single_particle || (s < THETA * r);
        if (accept_approx) {
            double denom = r + EPSILON;
            double f = G * mass * node->mass / (denom * denom * denom);
            local_fx += f * (-dx);
            local_fy += f * (-dy);
        } else {
            for (int i = 0; i < 4; i++) {
                if (node->child[i])
                    stack[sp++] = node->child[i];
            }
        }
    }

    *res_fx = local_fx;
    *res_fy = local_fy;
}
