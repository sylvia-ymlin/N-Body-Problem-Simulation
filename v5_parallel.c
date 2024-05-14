#include "core/domain.h"
#include "core/ds.h"
#include "core/morton.h"
#include "core/time_utils.h"
#include "simulation.h"
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double COINCIDENT_EPS = 1e-9;
static NodeArena arena = {NULL, 0, 0};
static int arena_N = 0;
static const double DOMAIN_PADDING_FRAC = 0.05;

typedef enum {
    SCHED_DYNAMIC = 0,
    SCHED_STATIC = 1,
    SCHED_GUIDED = 2,
} ScheduleKind;

static ScheduleKind load_schedule_kind(void);
static int load_schedule_chunk(void);
static const char* schedule_name(ScheduleKind kind);
static inline int is_leaf(TNode* node);
static inline int quadrant(double px, double py, double mx, double my);
static TNode* create_node(double LB, double RB, double DB, double UB,
                          NodeArena* a);
static void insert(TNode* node, int idx, ParticleSystem* sys, NodeArena* a);
static void compute_force_single(double pos_x, double pos_y, double mass,
                                 int particle_idx, TNode* root, double* res_fx,
                                 double* res_fy, double G, double THETA);

void compute_force_v5_parallel(ParticleSystem* sys, KernelConfig* config) {
    int N = sys->N;
    ScheduleKind schedule_kind = load_schedule_kind();
    int schedule_chunk = load_schedule_chunk();

    const double* x = sys->pos_x;
    const double* y = sys->pos_y;
    const double* m = sys->mass;
    double* fx_out = sys->fx;
    double* fy_out = sys->fy;

    static int domain_initialized = 0;
    static int domain_N = 0;
    static double x_min = 0.0, x_max = 0.0, y_min = 0.0, y_max = 0.0;
    if (!domain_initialized || domain_N != N) {
        domain_square_init_from_particles(x, y, N, DOMAIN_PADDING_FRAC, &x_min,
                                          &x_max, &y_min, &y_max);
        domain_initialized = 1;
        domain_N = N;
    } else {
        domain_square_expand_if_needed(x, y, N, DOMAIN_PADDING_FRAC, &x_min,
                                       &x_max, &y_min, &y_max);
    }

    z_order_sort(sys, x_min, x_max, y_min, y_max);

    double t_tree_start = sim_time_now();
    if (arena.buffer == NULL || arena_N != N) {
        if (arena.buffer)
            free_arena(&arena);
        init_arena(&arena, N * 100);
        arena_N = N;
    }
    reset_arena(&arena);
    TNode* root = create_node(x_min, x_max, y_min, y_max, &arena);
    for (int i = 0; i < N; i++)
        insert(root, i, sys, &arena);
    double t_tree_end = sim_time_now();

    double t_force_start = sim_time_now();
    double G = G_FACTOR / N;
#ifdef _OPENMP
    omp_sched_t omp_schedule = omp_sched_dynamic;
    if (schedule_kind == SCHED_STATIC)
        omp_schedule = omp_sched_static;
    else if (schedule_kind == SCHED_GUIDED)
        omp_schedule = omp_sched_guided;
    omp_set_schedule(omp_schedule, schedule_chunk);
#pragma omp parallel for schedule(runtime) num_threads(config->n_threads)
#endif
    for (int i = 0; i < N; i++) {
        compute_force_single(x[i], y[i], m[i], i, root, &fx_out[i], &fy_out[i],
                             G, config->theta_max);
    }
    double t_force_end = sim_time_now();

    printf("\n[HP] Schedule: %s,%d | Tree Build: %.4f s | Force Compute: %.4f s\n",
           schedule_name(schedule_kind), schedule_chunk,
           t_tree_end - t_tree_start, t_force_end - t_force_start);
}

static ScheduleKind load_schedule_kind(void) {
    const char* raw = getenv("NBODY_OMP_SCHEDULE");
    if (!raw || raw[0] == '\0')
        return SCHED_DYNAMIC;

    char buf[32];
    size_t n = strlen(raw);
    if (n >= sizeof(buf))
        n = sizeof(buf) - 1;
    for (size_t i = 0; i < n; i++)
        buf[i] = (char)tolower((unsigned char)raw[i]);
    buf[n] = '\0';

    if (strcmp(buf, "static") == 0)
        return SCHED_STATIC;
    if (strcmp(buf, "guided") == 0)
        return SCHED_GUIDED;
    return SCHED_DYNAMIC;
}

static int load_schedule_chunk(void) {
    const char* raw = getenv("NBODY_OMP_CHUNK");
    if (!raw || raw[0] == '\0')
        return 128;
    int chunk = atoi(raw);
    return chunk > 0 ? chunk : 128;
}

static const char* schedule_name(ScheduleKind kind) {
    switch (kind) {
    case SCHED_STATIC:
        return "static";
    case SCHED_GUIDED:
        return "guided";
    default:
        return "dynamic";
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
    double px = sys->pos_x[idx], py = sys->pos_y[idx], m = sys->mass[idx];

    if (node->particle_idx == -1 && is_leaf(node) && node->mass == 0) {
        node->particle_idx = idx;
        node->mass = m;
        node->pos_x = px;
        node->pos_y = py;
        return;
    }

    if (is_leaf(node)) {
        int old = node->particle_idx;
        if (old != -1) {
            if (fabs(px - sys->pos_x[old]) < COINCIDENT_EPS &&
                fabs(py - sys->pos_y[old]) < COINCIDENT_EPS) {
                double mt = node->mass + m;
                node->pos_x = (node->pos_x * node->mass + px * m) / mt;
                node->pos_y = (node->pos_y * node->mass + py * m) / mt;
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

    double mt = node->mass + m;
    node->pos_x = (node->pos_x * node->mass + px * m) / mt;
    node->pos_y = (node->pos_y * node->mass + py * m) / mt;
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

        if (node->particle_idx != -1 || (s < THETA * r)) {
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
