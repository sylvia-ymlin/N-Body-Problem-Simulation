#include "morton.h"

// Magic numbers for spreading 32-bit integer into 64-bit integer
// leaving zeros in even positions.
// 0000...0000 dcba -> 0d0c0b0a
static inline uint64_t split_by_2(unsigned int a) {
  uint64_t x = a & 0xffffffff;
  x = (x | x << 16) & 0x0000ffff0000ffff;
  x = (x | x << 8) & 0x00ff00ff00ff00ff;
  x = (x | x << 4) & 0x0f0f0f0f0f0f0f0f;
  x = (x | x << 2) & 0x3333333333333333;
  x = (x | x << 1) & 0x5555555555555555;
  return x;
}

static inline uint64_t morton_encode_magicbits(unsigned int x, unsigned int y) {
  return (split_by_2(x) | (split_by_2(y) << 1));
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

void z_order_sort(double *pos_x, double *pos_y, double *mass, double *vx,
                  double *vy, double *brightness, int N, double LB, double RB,
                  double DB, double UB) {

  SortEntry *entries = (SortEntry *)malloc(N * sizeof(SortEntry));
  if (!entries)
    return;

  // Normalization factors
  double width = RB - LB;
  double height = UB - DB;
  // Map to 32-bit integer (max for split_by_2)
  double scale_x = (double)((1ULL << 32) - 1) / width;
  double scale_y = (double)((1ULL << 32) - 1) / height;

  // 1. Compute codes
  for (int i = 0; i < N; i++) {
    unsigned int ix = (unsigned int)((pos_x[i] - LB) * scale_x);
    unsigned int iy = (unsigned int)((pos_y[i] - DB) * scale_y);
    entries[i].index = i;
    entries[i].code = morton_encode_magicbits(ix, iy);
  }

  // 2. Sort indices based on code
  qsort(entries, N, sizeof(SortEntry), compare_entries);

  // 3. Permute arrays
  double *temp = (double *)malloc(N * sizeof(double));
  if (!temp) {
    free(entries);
    return;
  }

// Helper macro for permutation
#define PERMUTE(arr)                                                           \
  for (int i = 0; i < N; i++)                                                  \
    temp[i] = arr[entries[i].index];                                           \
  memcpy(arr, temp, N * sizeof(double));

  PERMUTE(pos_x);
  PERMUTE(pos_y);
  PERMUTE(mass);
  PERMUTE(vx);
  PERMUTE(vy);
  PERMUTE(brightness);

#undef PERMUTE

  free(temp);
  free(entries);
}
