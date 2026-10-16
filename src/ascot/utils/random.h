/**
 * @file random.h
 * Random number generator interface.
 *
 * This module uses random123 library which provides counter based random number
 * generators (as opposed to traditional state based). This means that every
 * sampled random number is fully defined by the seed, marker ID and the RNG
 * counter index (the counter is marker specific and advanced by this module).
 * Therefore each marker always gets the same sequence of random numbers when
 * using the same seed (and marker ID) irrespective of simulation size or any
 * race conditions. This makes the simulation reproducible and deterministic
 * as long as there are no other sources of randomness.
 */
#ifndef RANDOM_H
#define RANDOM_H

#include "ascot.h"
#include "random123/philox.h"
#include <stdlib.h>

/**
 * Initialize random number generator.
 *
 * This should only be called between simulations.
 *
 * @param seed Seed for random number generator.
 */
void random_seed(uint64_t seed);

/**
 * Sample random numbers from uniform distribution (0.0, 1.0].
 *
 * This uses Philox4x32 random number generator which creates 32bit random
 * numbers in blocks of four.
 *
 * @param rnd Array of generated random numbers.
 * @param nrnd Number of random numbers to generate.
 * @param marker_id Array of marker IDs.
 * @param rng_count Array of random number counters. This value is advanced as
 *        many times as needed based on the number of random numbers requested
 *        and underlying algorithm.
 * @param index Index of the current marker.
 */
void random_uniform(
    real *rnd, size_t nrnd, const size_t *marker_id, size_t *rng_count,
    size_t index);

/**
 * Sample random numbers from Normal distribution.
 *
 * This uses Philox4x32 random number generator which creates 32bit random
 * numbers in blocks of four.
 *
 * @param rnd Array of generated random numbers.
 * @param nrnd Number of random numbers to generate.
 * @param marker_id Array of marker IDs.
 * @param rng_count Array of random number counters. This value is advanced as
 *        many times as needed based on the number of random numbers requested
 *        and underlying algorithm.
 * @param index Index of the current marker.
 */
void random_norm(
    real *rnd, size_t nrnd, const size_t *marker_id, size_t *rng_count,
    size_t index);

#endif
