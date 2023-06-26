/**
 * @file libascot_mem.h
 * @brief Header file for libascot_mem.c
 */
#ifndef LIBASCOT_MEM_H_
#define LIBASCOT_MEM_H_

#include "particle.h"

input_particle* libascot_allocate_input_particles(int nmrk);
particle_state* libascot_allocate_particle_states(int nmrk);
real* libascot_allocate_reals(size_t size);

void libascot_deallocate(void *arr);

#endif
