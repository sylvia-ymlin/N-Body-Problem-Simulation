#ifndef MORTON_H
#define MORTON_H

#include "types.h"

// Reorder particles by Morton code within the given bounding box.
void z_order_sort(ParticleSystem* sys, double LB, double RB, double DB,
                  double UB);

#endif
