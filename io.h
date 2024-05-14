#ifndef IO_H
#define IO_H

#include "types.h"

ParticleSystem io_read_particles(const char* filename, int N);
void io_write_result(const char* filename, ParticleSystem* sys);
void io_free_particles(ParticleSystem* sys);

#endif
