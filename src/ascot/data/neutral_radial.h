/**
 * @file neutral_radial.h
 * Radial neutral data with linear interpolation.
 */
#ifndef NEUTRAL_RADIAL_H
#define NEUTRAL_RADIAL_H

#include "defines.h"
#include "neutral.h"
#include "parallel.h"

/**
 * Initialize the linearly interpolated radial neutral data.
 *
 * @param neutral The struct to initialize.
 * @param n Number of neutral species.
 * @param nrho Number of rho grid points in the data.
 * @param rholim Limits of the uniform rho grid [1].
 * @param density Density of the neutral species (in same order as they are
 *        listed in the plasma input) [m^-3].
 *        Layout is (neutral, rhoi) = [neutral*nrho + i] (C order).
 * @param temperature Temperature of the neutral species (in same order as they
 *        are listed in the plasma input) [J].
 *        Layout is (neutral, rhoi) = [neutral*nrho + i] (C order).
 *
 * @return Zero if initialization succeeded.
 */
int NeutralRadial_init(
    NeutralRadial *neutral, size_t n, size_t nrho, real rholim[2],
    real density[n * nrho], real temperature[n * nrho]);

/**
 * Free allocated resources.
 *
 * @param neutral The struct whose fields are deallocated.
 */
void NeutralRadial_free(NeutralRadial *neutral);

/**
 * Offload data to the accelerator.
 *
 * @param neutral The struct to offload.
 */
void NeutralRadial_offload(NeutralRadial *neutral);

DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate density of each neutral species.
 *
 * The evaluated values are in the same order as ions in the plasma data.
 *
 * @param density Evaluated density [m^-3].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t NeutralRadial_eval_density(
    real *density, real rho, NeutralRadial *neutral);

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
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t NeutralRadial_eval_temperature(
    real *temperature, real rho, NeutralRadial *neutral);

#endif
