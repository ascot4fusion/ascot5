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

#include "atomic_reactions.h"
#include "data/atomic.h"
#include "data/marker.h"
#include "data/neutral.h"
#include "data/plasma.h"
#include "defines.h"
#include "utils/random.h"

/**
 * Check and apply atomic reactions.
 *
 * @param mrk Markers to which atomic reactions are applied.
 * @param h Time step [s].
 * @param plasma Plasma data.
 * @param neutral Neutral data.
 * @param atomic Atomic reaction data.
 */
void atomic_reactions_eval_go(
    MarkerGyroOrbit *mrk, real *h, Plasma *plasma, Neutral *neutral,
    Atomic *atomic);

GPU_DECLARE_TARGET_SIMD_UNIFORM(atomic)
/**
 * Determines atomic reaction rates
 *
 * The reaction rates are determined based on the charge state. The terms
 * "ionization" and "recombination" are used loosely, where "ionization" refers
 * to all charge-increasing reactions, and "recombination" to all
 * charge-decreasing reactions.
 *
 * @param rate_ionization Evaluated reaction rate for ionization [1/s].
 * @param rate_recombination Evaluated reaction rate for recombination [1/s].
 * @param z_1 Test particle atomic number.
 * @param a_1 Test particle mass number.
 * @param m_1 Test particle mass [kg].
 * @param z_2 Atomic numbers of background plasma species.
 * @param a_2 Mass numbers of background plasma species.
 * @param atomic Atomic reaction data.
 * @param q Test particle charge state.
 * @param E Test particle energy [J].
 * @param N_pls_spec Number of background plasma species.
 * @param N_ntl_spec Number of background neutral species.
 * @param T Background plasma temperature [J].
 * @param T_0 Background neutral temperature [J].
 * @param n Density of background plasma species [m^-3].
 * @param n_0 Density of background neutrals [m^-3].
 *
 * @return zero if evaluation succeeded
 */
err_t atomic_reactions_rates(
    real *rate_ionization, real *rate_recombination, int z_1, int a_1, real m_1,
    const int *z_2, const int *a_2, Atomic *atomic, int q,
    real E, int N_pls_spec, int N_ntl_spec, real *T, real *T_0, real *n,
    real *n_0);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD
/**
 * Apply atomic reactions.
 *
 * Change (or keep intact) marker charge state based on the reaction rates.
 *
 * @param charge Pointer to the test particle charge which is modified if a
 *        reaction occurs.
 * @param znum Test particle atomic number.
 * @param h Time step [s].
 * @param rate_ionization Combined reaction rate for all ionization processes
 *        [1/s].
 * @param rate_recombination Combined reaction rate for all recombination
 *        processes [1/s].
 * @param rnd Uniformly (0, 1] distributed random numbers.
 */
void atomic_reactions_apply(
    int *charge, int znum, real h, real rate_ionization,
    real rate_recombination, real rnd[2]);
GPU_DECLARE_TARGET_SIMD_END

#endif
