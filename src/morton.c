#include "morton.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Simple bit interleaving for clarity
static inline uint64_t morton_encode(unsigned int x, unsigned int y) {
  uint64_t answer = 0;
  for (uint64_t i = 0; i < 32; i++) { // Interleave bits
    // Use x[i] -> answer[2*i]
    // Use y[i] -> answer[2*i+1]
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

  // Normalization factors
  double width = RB - LB;
  double height = UB - DB;
  // Map to 32-bit integer (max for split_by_2)
  double scale_x = (double)((1ULL << 32) - 1) / width;
  const double scale_y = (double)((1ULL << 32) - 1) / height;

  // 1. Compute codes
  for (int i = 0; i < N; i++) {
    unsigned int ix = (unsigned int)((sys->pos_x[i] - LB) * scale_x);
    unsigned int iy = (unsigned int)((sys->pos_y[i] - DB) * scale_y);
    entries[i].index = i;
    entries[i].code = morton_encode(ix, iy);
  }

  // 2. Sort indices based on code
  qsort(entries, N, sizeof(SortEntry), compare_entries);

  // 3. Permute arrays
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
  // PERMUTE(sys->fx); // No need to permute forces if they are recalculated
  // every step, but safer to do so or zero them? Standard practice: forces are
  // scratchpad, but if we permute mid-step, we might need them. However,
  // sorting usually happens at start or between steps. Given user request
  // "force only depends on pos", per-step sort implies we might want to permute
  // fx/fy if we do it in middle of step. But usually we sort before force calc.
  // I will permute fx/fy to be safe/consistent, though usually they are reset.
  PERMUTE(sys->fx);
  PERMUTE(sys->fy);

#undef PERMUTE

  free(temp);
  free(entries);
}
