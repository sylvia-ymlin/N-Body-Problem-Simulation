#ifndef TIME_UTILS_H
#define TIME_UTILS_H

#ifdef _OPENMP
#include <omp.h>
#else
#include <sys/time.h>
#endif

static inline double sim_time_now() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
#endif
}

#endif
