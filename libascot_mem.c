#include "libascot_mem.h"
#include <stdlib.h>

/**
 *  @brief A routine to allocate an array of particles
 *
 */
input_particle* libascot_allocate_input_particles(int number_of_particles){
	return (input_particle*) malloc(number_of_particles * sizeof(input_particle) );
}

/**
 *  @brief A routine to allocate an array of reals
 *
 */
real* libascot_allocate_reals(size_t size){
  return ( real* ) malloc( size * sizeof(real) );
}

/**
 *  @brief A wrapper to C stdlib free()
 *
 */
void libascot_deallocate(void *arr){
  free( arr );
}
