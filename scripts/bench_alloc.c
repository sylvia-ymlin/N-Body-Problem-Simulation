#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#define N_NODES 100000
#define N_ITERATIONS 100

typedef struct TNode {
  double data[4];
  struct TNode *child[4];
} TNode;

double get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

// --- Naive Approach ---
TNode **naive_ptrs;
void bench_naive() {
  for (int iter = 0; iter < N_ITERATIONS; iter++) {
    // Allocate
    for (int i = 0; i < N_NODES; i++) {
      naive_ptrs[i] = (TNode *)malloc(sizeof(TNode));
    }
    // Free
    for (int i = 0; i < N_NODES; i++) {
      free(naive_ptrs[i]);
    }
  }
}

// --- Arena Approach ---
typedef struct Arena {
  TNode *buffer;
  size_t used;
} Arena;
Arena my_arena;

void bench_arena() {
  for (int iter = 0; iter < N_ITERATIONS; iter++) {
    // Allocate (Reset)
    my_arena.used = 0;

    // Simulate construction
    for (int i = 0; i < N_NODES; i++) {
      TNode *node = &my_arena.buffer[my_arena.used++];
      // Optional: simulate touch to ensure page fault/cache effect is fair?
      // node->data[0] = 0;
    }
  }
}

int main() {
  naive_ptrs = malloc(N_NODES * sizeof(TNode *));
  my_arena.buffer = malloc(N_NODES * sizeof(TNode));

  // Warmup
  bench_naive();
  bench_arena();

  // Bench Naive
  double start = get_time();
  bench_naive();
  double end = get_time();
  printf("Naive Malloc/Free: %.4f seconds\n", end - start);

  // Bench Arena
  start = get_time();
  bench_arena();
  end = get_time();
  printf("Linear Arena:      %.4f seconds\n", end - start);

  return 0;
}
