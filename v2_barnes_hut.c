#include "core/domain.h"
#include "core/ds.h"
#include "simulation.h"
#include <math.h>
#include <stdlib.h>

static const double COINCIDENT_EPS = 1e-9;
static const double DOMAIN_PADDING_FRAC = 0.05;

static inline int is_leaf(TNode* node);
static inline int quadrant(double px, double py, double mx, double my);
static TNode* create_node(double LB, double RB, double DB, double UB);
static void free_tree(TNode* node);
static void insert(TNode* node, int idx, ParticleSystem* sys);
static TNode* build_tree(ParticleSystem* sys, int N, double x_min, double x_max,
                         double y_min, double y_max);
static void accumulate_force_rec(TNode* node, int idx, ParticleSystem* sys,
                                 double* fx, double* fy, double theta,
                                 double G);

void compute_force_v2_barnes_hut(ParticleSystem* sys, KernelConfig* config) {
    const int N = sys->N;
    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* mass_arr = sys->mass;
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

    const double G = G_FACTOR / N;
    const double theta = config->theta_max;

    TNode* root = build_tree(sys, N, domain.x_min, domain.x_max, domain.y_min,
                             domain.y_max);
    for (int i = 0; i < N; i++) {
        fx_out[i] = 0;
        fy_out[i] = 0;
        accumulate_force_rec(root, i, sys, &fx_out[i], &fy_out[i], theta, G);
    }

    free_tree(root);
}

static TNode* build_tree(ParticleSystem* sys, int N, double x_min, double x_max,
                         double y_min, double y_max) {
    TNode* root = create_node(x_min, x_max, y_min, y_max);
    for (int i = 0; i < N; i++)
        insert(root, i, sys);
    return root;
}

static void accumulate_force_rec(TNode* node, int idx, ParticleSystem* sys,
                                 double* fx, double* fy, double theta,
                                 double G) {
    if (!node || node->mass == 0 || node->particle_idx == idx)
        return;

    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* mass_arr = sys->mass;

    double dx = node->pos_x - x[idx];
    double dy = node->pos_y - y[idx];
    double r = sqrt(dx * dx + dy * dy);
    double cell_width = node->x_max - node->x_min;

    const int accept_approx = is_leaf(node) || (cell_width / r < theta);
    if (accept_approx) {
        double denom = r + EPSILON;
        double f = G * mass_arr[idx] * node->mass / (denom * denom * denom);
        *fx += f * dx;
        *fy += f * dy;
    } else {
        for (int i = 0; i < 4; i++)
            accumulate_force_rec(node->child[i], idx, sys, fx, fy, theta, G);
    }
}

static void insert(TNode* node, int idx, ParticleSystem* sys) {
    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* mass_arr = sys->mass;
    double px = x[idx], py = y[idx], mass = mass_arr[idx];

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
                (oq & 2) ? my : node->y_min, (oq & 2) ? node->y_max : my);
            insert(node->child[oq], old, sys);

            node->mass = old_m;
            node->pos_x = old_px;
            node->pos_y = old_py;
        }
    }

    double mx = (node->x_min + node->x_max) * 0.5,
           my = (node->y_min + node->y_max) * 0.5;
    int q = quadrant(px, py, mx, my);
    if (!node->child[q]) {
        node->child[q] =
            create_node((q & 1) ? mx : node->x_min, (q & 1) ? node->x_max : mx,
                        (q & 2) ? my : node->y_min, (q & 2) ? node->y_max : my);
    }
    insert(node->child[q], idx, sys);

    double mt = node->mass + mass;
    node->pos_x = (node->pos_x * node->mass + px * mass) / mt;
    node->pos_y = (node->pos_y * node->mass + py * mass) / mt;
    node->mass = mt;
}

static TNode* create_node(double LB, double RB, double DB, double UB) {
    TNode* node = (TNode*)malloc(sizeof(TNode));
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

static void free_tree(TNode* node) {
    if (!node)
        return;
    for (int i = 0; i < 4; i++)
        free_tree(node->child[i]);
    free(node);
}

static inline int is_leaf(TNode* node) {
    return node->child[0] == NULL && node->child[1] == NULL &&
           node->child[2] == NULL && node->child[3] == NULL;
}

static inline int quadrant(double px, double py, double mx, double my) {
    return (px > mx) + 2 * (py > my);
}
