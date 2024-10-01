/**
 * @file bbnbi5.h
 * @brief Functions to execute bbnbi externally
 */
#ifndef BBNBI5_H
#define BBNBI5_H

#include "ascot5.h"
#include "simulate.h"

void bbnbi_simulate(
    sim_data* sim, int nprt, real t1, real t2, particle_state** p);

#endif
