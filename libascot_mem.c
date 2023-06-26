/**
 * @file libascot_mem.c
 * @brief Provides memory de/allocation routines for the libascot
 */
#include <stdlib.h>
#include "libascot_mem.h"
#include "particle.h"

/**
 * @brief A routine to allocate an array of input particles
 *
 * @param nmrk number of markers for which space is allocated
 *
 * @return pointer to the allocated array
 */
input_particle* libascot_allocate_input_particles(int nmrk) {
    return (input_particle*) malloc(nmrk * sizeof(input_particle) );
}

/**
 * @brief A routine to allocate an array of particle states
 *
 * @param nmrk number of markers for which space is allocated
 *
 * @return pointer to the allocated array
 */
particle_state* libascot_allocate_particle_states(int nmrk) {
    return (particle_state*) malloc(nmrk * sizeof(particle_state) );
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
