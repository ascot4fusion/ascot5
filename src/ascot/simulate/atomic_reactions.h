/**
 * @file atomic_reactions.h
 * Atomic reactions model
 *
 * The atomic module models atomic reactions for the simulated test particles.
 * Using atomic reaction data supplied by the asigma (atomicsigma) helper
 * module, reaction probabilities are calculated and tested against a random
 * number. If a reaction occurs, the particle charge state is changed.
 *
 * The terms ionization and recombination are used loosely in variable names.
 * Ionization (ion) stands for all charge-increasing reactions, and
 * recombination (rec) for all charge-decreasing reactions.
 */
#ifndef ATOMIC_REACTIONS_H
#define ATOMIC_REACTIONS_H

#include "defines.h"
#include "atomic_reactions.h"
#include "atomic.h"
#include "neutral.h"
#include "particle.h"
#include "plasma.h"
#include "random.h"

/**
 * @brief Determine if atomic reactions occur during time-step and change charge
 *
 * @param p fo struct
 * @param h time-steps from NSIMD markers
 * @param p_data pointer to plasma data
 * @param n_data pointer to neutral data
 * @param r_data pointer to random-generator data
 * @param asigmadata pointer to atomic reaction data
 */
#ifndef GPU
#pragma omp declare target
#endif
void atomic_fo(
    particle_simd_fo *p, real *h, Plasma *plasma, Neutral *neutral,
    random_data *r_data, Atomic *atomic);
#ifndef GPU
#pragma omp end declare target
#endif
#endif
