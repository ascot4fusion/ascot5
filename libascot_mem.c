#include "libascot_mem.h"
#include <stdlib.h>

/**
 *  @brief A routine to allocate an array of particles
 *
 */
input_particle* libascot_allocate_input_particles(int number_of_particles){
	return (input_particle*) malloc(number_of_particles * sizeof(input_particle) );
}
