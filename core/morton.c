#include "morton.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Interleave x/y bits into a Morton code.
static inline uint64_t morton_encode(unsigned int x, unsigned int y) {
  uint64_t answer = 0;
  for (uint64_t i = 0; i < 32; i++) {
    uint64_t x_bit = (x >> i) & 1;
    uint64_t y_bit = (y >> i) & 1;
    answer |= (x_bit << (2 * i)) | (y_bit << (2 * i + 1));
  }
  return answer;
}

typedef struct SortEntry {
  int index;
  uint64_t code;
} SortEntry;

int compare_entries(const void *a, const void *b) {
  uint64_t ca = ((SortEntry *)a)->code;
  uint64_t cb = ((SortEntry *)b)->code;
  if (ca < cb)
    return -1;
  if (ca > cb)
    return 1;
  return 0;
}

void z_order_sort(ParticleSystem *sys, double LB, double RB, double DB,
                  double UB) {
  int N = sys->N;
  SortEntry *entries = (SortEntry *)malloc(N * sizeof(SortEntry));
  if (!entries)
    return;

  double width = RB - LB;
  double height = UB - DB;
  // Normalize coordinates into a 32-bit integer grid before bit interleaving.
  double scale_x = (double)((1ULL << 32) - 1) / width;
  const double scale_y = (double)((1ULL << 32) - 1) / height;

  for (int i = 0; i < N; i++) {
    unsigned int ix = (unsigned int)((sys->pos_x[i] - LB) * scale_x);
    unsigned int iy = (unsigned int)((sys->pos_y[i] - DB) * scale_y);
    entries[i].index = i;
    entries[i].code = morton_encode(ix, iy);
  }

  qsort(entries, N, sizeof(SortEntry), compare_entries);

  double *temp = (double *)malloc(N * sizeof(double));
  if (!temp) {
    free(entries);
    return;
  }

#define PERMUTE(arr)                                                           \
  for (int i = 0; i < N; i++)                                                  \
    temp[i] = arr[entries[i].index];                                           \
  memcpy(arr, temp, N * sizeof(double));

  PERMUTE(sys->pos_x);
  PERMUTE(sys->pos_y);
  PERMUTE(sys->mass);
  PERMUTE(sys->vx);
  PERMUTE(sys->vy);
  // Force arrays are scratch state, but keeping them aligned avoids subtle
  // mismatches if sorting happens between force evaluations.
  PERMUTE(sys->fx);
  PERMUTE(sys->fy);

#undef PERMUTE

  free(temp);
  free(entries);
}
