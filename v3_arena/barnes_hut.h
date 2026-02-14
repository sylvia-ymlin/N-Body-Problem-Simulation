#ifndef BARNES_HUT_H
#define BARNES_HUT_H

#include "ds.h"

void barnes_hut(double *px, double *py, double *mass, int N, double *fx,
                double *fy, double THETA_MAX, NodeArena *arena);

#endif
