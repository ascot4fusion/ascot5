/**
 * @file libascot_mem.h
 * @brief Header file for libascot_mem.c
 *
 *  Created on: Nov 25, 2021
 *      Author: sjjamsa
 */
#ifndef LIBASCOT_MEM_H_
#define LIBASCOT_MEM_H_

#include "particle.h"

input_particle* libascot_allocate_input_particles(int number_of_particles);

real* libascot_allocate_reals(size_t size);

void libascot_deallocate(void *arr);

#endif
