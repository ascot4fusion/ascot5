/*
 * libascot_mem.h
 *
 *  Created on: Nov 25, 2021
 *      Author: sjjamsa
 */

#ifndef LIBASCOT_MEM_H_
#define LIBASCOT_MEM_H_

#include "particle.h"

/**
 *  @brief A routine to allocate an array of particles
 *
 */
input_particle* libascot_allocate_input_particles(int number_of_particles);


#endif /* LIBASCOT_MEM_H_ */
