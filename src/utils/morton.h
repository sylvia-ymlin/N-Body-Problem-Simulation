#ifndef MORTON_H
#define MORTON_H

#include "../core/types.h"

/**
 * Z-Order Sorting (Morton Code Sorting)
 * Reorders particles in the system based on their spatial position.
 *
 * @param sys: Pointer to ParticleSystem
 * @param LB, RB, DB, UB: Bounding box (Left, Right, Down, Up)
 */
void z_order_sort(ParticleSystem *sys, double LB, double RB, double DB,
                  double UB);

#endif
