#ifndef MORTON_H
#define MORTON_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Z-Order Sorting (Morton Code Sorting)
 * Reorders all particle arrays based on their spatial position to improve cache
 * locality.
 *
 * @param pos_x, pos_y, mass, vx, vy, brightness: Particle data arrays to be
 * permuted
 * @param N: Number of particles
 * @param LB, RB, DB, UB: Bounding box of the simulation region
 */
void z_order_sort(double *pos_x, double *pos_y, double *mass, double *vx,
                  double *vy, int N, double LB, double RB, double DB,
                  double UB);

#endif
