/**
 * @file neutral.h
 * Neutral interface.
 *
 * Provides functions and datatypes for getting neutral species information and
 * for interpolating background neutral density and temperature.
 *
 * The neutral data is closely related to the plasma data because it is assumed
 * that the species in the neutral data are the same, and in same order, as in
 * the plasma data.
 */
#ifndef NEUTRAL_H
#define NEUTRAL_H

#include "defines.h"
#include "neutral_data.h"
#include "parallel.h"

/**
 * Free allocated resources.
 *
 * @param neutral The struct whose fields are deallocated.
 */
void Neutral_free(Neutral *neutral);

/**
 * Offload data to the accelerator.
 *
 * @param neutral The struct to offload.
 */
void Neutral_offload(Neutral *neutral);

DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate density of each neutral species.
 *
 * Depending on the implementation, either the normalized poloidal flux
 * coordinate or the cylindrical coordinate is used in the interpolation.
 *
 * The evaluated values are in the same order as ions in the plasma data.
 *
 * @param density Evaluated density [m^-3].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Neutral_eval_density(
    real *density, real rho, real r, real phi, real z, real t,
    Neutral *neutral);

DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate temperature of each neutral species.
 *
 * Depending on the implementation, either the normalized poloidal flux
 * coordinate or the cylindrical coordinate is used in the interpolation.
 *
 * The evaluated values are in the same order as ions in the plasma data.
 *
 * @param temperature Evaluated temperature [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Neutral_eval_temperature(
    real *temperature, real rho, real r, real phi, real z, real t,
    Neutral *neutral);

DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Get the number of neutral species.
 *
 * @param neutral The neutral data.
 *
 * @return The number of neutral species.
 */
size_t Neutral_get_n_species(Neutral *neutral);

#endif
