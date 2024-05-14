#include "morton.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/** Interleave x and y bits into a 64-bit Morton (Z-order) code.
 * ----------------------------------------------------------------- */
static uint64_t morton_encode(unsigned int x, unsigned int y) {
    uint64_t code = 0;
    for (uint64_t i = 0; i < 32; i++) {
        uint64_t x_bit = (x >> i) & 1;
        uint64_t y_bit = (y >> i) & 1;
        code |= (x_bit << (2 * i)) | (y_bit << (2 * i + 1));
    }
    return code;
}

typedef struct {
    int      index;
    uint64_t code;
} SortEntry;

static int compare_entries(const void* a, const void* b) {
    uint64_t ca = ((SortEntry*)a)->code;
    uint64_t cb = ((SortEntry*)b)->code;
    if (ca < cb) return -1;
    if (ca > cb) return  1;
    return 0;
}

/** Reorder all particle arrays by Z-order curve within the given bounding box.
 * Nearby particles in space end up nearby in memory after this sort.
 * ----------------------------------------------------------------- */
void z_order_sort(ParticleSystem* sys, double LB, double RB, double DB, double UB) {
    int N = sys->N;
    SortEntry* entries = (SortEntry*)malloc(N * sizeof(SortEntry));
    if (!entries) return;

    double scale_x = (double)((1ULL << 32) - 1) / (RB - LB);
    double scale_y = (double)((1ULL << 32) - 1) / (UB - DB);

    for (int i = 0; i < N; i++) {
        unsigned int ix = (unsigned int)((sys->pos_x[i] - LB) * scale_x);
        unsigned int iy = (unsigned int)((sys->pos_y[i] - DB) * scale_y);
        entries[i].index = i;
        entries[i].code  = morton_encode(ix, iy);
    }

    qsort(entries, N, sizeof(SortEntry), compare_entries);

    /* Permute each particle array into the new order */
    double* temp = (double*)malloc(N * sizeof(double));
    if (!temp) { free(entries); return; }

    for (int i = 0; i < N; i++) temp[i] = sys->pos_x[entries[i].index];
    memcpy(sys->pos_x, temp, N * sizeof(double));

    for (int i = 0; i < N; i++) temp[i] = sys->pos_y[entries[i].index];
    memcpy(sys->pos_y, temp, N * sizeof(double));

    for (int i = 0; i < N; i++) temp[i] = sys->mass[entries[i].index];
    memcpy(sys->mass, temp, N * sizeof(double));

    for (int i = 0; i < N; i++) temp[i] = sys->vx[entries[i].index];
    memcpy(sys->vx, temp, N * sizeof(double));

    for (int i = 0; i < N; i++) temp[i] = sys->vy[entries[i].index];
    memcpy(sys->vy, temp, N * sizeof(double));

    for (int i = 0; i < N; i++) temp[i] = sys->fx[entries[i].index];
    memcpy(sys->fx, temp, N * sizeof(double));

    for (int i = 0; i < N; i++) temp[i] = sys->fy[entries[i].index];
    memcpy(sys->fy, temp, N * sizeof(double));

    free(temp);
    free(entries);
}
