/**
 * @file libascot_mem.c
 * @brief Provides memory de/allocation routines for the libascot
 */
#include "libascot_mem.h"
#include <stdlib.h>


/**
 * @brief A routine to allocate an array of particles
 *
 * @param number_of_particles number of markers for which space is allocated
 *
 * @return pointer to the allocated array
 */
input_particle* libascot_allocate_input_particles(int number_of_particles) {
    return (input_particle*) malloc(number_of_particles * sizeof(input_particle) );
}


/**
 * @brief A routine to allocate an array of reals
 *
 * @param size size of the array
 *
 * @return pointer to the allocated array
 */
real* libascot_allocate_reals(size_t size) {
    return ( real* ) malloc( size * sizeof(real) );
}


/**
 * @brief A wrapper to C stdlib free()
 *
 * @param arr array to be freed
 */
void libascot_deallocate(void *arr) {
    free(arr);
}
