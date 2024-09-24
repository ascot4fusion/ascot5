/**
 * @file bbnbi5.h
 * @brief Functions to execute bbnbi externally
 */
#ifndef BBNBI5_H
#define BBNBI5_H

#include "ascot5.h"
#include "simulate.h"

void bbnbi_simulate(
    sim_offload_data* sim, int nprt, real t1, real t2, real* B_offload_array,
    real* plasma_offload_array, real* neutral_offload_array,
    real* wall_offload_array, int* wall_int_offload_array,
    real* asigma_offload_array, real* nbi_offload_array, particle_state** p,
    real* diag_offload_array);

#endif
