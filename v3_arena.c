#include "core/domain.h"
#include "core/ds.h"
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
static void compute_force_rec(TNode* node, int idx, ParticleSystem* sys,
                              double* fx, double* fy, double THETA, double G);

void compute_force_v3_arena(ParticleSystem* sys, KernelConfig* config) {
    if (arena.buffer == NULL || arena_N != sys->N) {
        if (arena.buffer)
            free_arena(&arena);
        init_arena(&arena, sys->N * ARENA_NODE_FACTOR);
        arena_N = sys->N;
    }
    reset_arena(&arena);

    int N = sys->N;
    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
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

    TNode* root = create_node(domain.x_min, domain.x_max, domain.y_min,
                              domain.y_max, &arena);
    for (int i = 0; i < N; i++)
        insert(root, i, sys, &arena);

    double G = G_FACTOR / N;
    double theta = config->theta_max;

    for (int i = 0; i < N; i++) {
        fx_out[i] = 0;
        fy_out[i] = 0;
        compute_force_rec(root, i, sys, &fx_out[i], &fy_out[i], theta, G);
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

static void compute_force_rec(TNode* node, int idx, ParticleSystem* sys,
                              double* fx, double* fy, double THETA, double G) {
    if (!node || node->mass == 0 || node->particle_idx == idx)
        return;

    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* m = sys->mass;

    double dx = node->pos_x - x[idx];
    double dy = node->pos_y - y[idx];
    double r = sqrt(dx * dx + dy * dy);
    double s = node->x_max - node->x_min;

    if (is_leaf(node) || (s / r < THETA)) {
        double denom = r + EPSILON;
        double f = G * m[idx] * node->mass / (denom * denom * denom);
        *fx += f * dx;
        *fy += f * dy;
    } else {
        for (int i = 0; i < 4; i++)
            compute_force_rec(node->child[i], idx, sys, fx, fy, THETA, G);
    }
}
