#ifndef IO_H
#define IO_H

#include "types.h"

// Reads particles from a file. Allocates memory in ParticleSystem.
ParticleSystem io_read_particles(const char *filename, int N);

// Appends current frame (pos, mass) to the movie file.
void io_write_frame(const char *filename, ParticleSystem *sys);

// Writes final state (pos, mass, vel) to the result file.
void io_write_result(const char *filename, ParticleSystem *sys);

// Frees memory allocated in ParticleSystem.
void io_free_particles(ParticleSystem *sys);

#endif
