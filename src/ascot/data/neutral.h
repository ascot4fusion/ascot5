/**
 * @file neutral.h
 * Neutral interface.
 *
 * Provides functions and datatypes for getting neutral species information and
 * for interpolating background neutral density and temperature.
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
 * Evaluate density of a neutral species.
 *
 * Depending on the implementation, either the normalized poloidal flux
 * coordinate or the cylindrical coordinate is used in the interpolation.
 *
 * @param density Evaluated density [m^-3].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param i_species Index of the requested species.
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Neutral_eval_density(
    real density[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Neutral *neutral);

DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate temperature of a neutral species.
 *
 * Depending on the implementation, either the normalized poloidal flux
 * coordinate or the cylindrical coordinate is used in the interpolation.
 *
 * @param temperature Evaluated temperature [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param i_species Index of the requested species.
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Neutral_eval_temperature(
    real temperature[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Neutral *neutral);

GPU_DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate density and temperature for all neutral species.
 *
 * Depending on the implementation, either the normalized poloidal flux
 * coordinate or the cylindrical coordinate is used in the interpolation.
 *
 * @param density Evaluated density [m^-3].
 * @param temperature Evaluated temperature [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r Radial coordinate of the query point [m].
 * @param phi Toroidal coordinate of the query point [rad].
 * @param z Axial coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Neutral_eval_nT(
    real *density, real *temperature, real rho, real r, real phi, real z,
    real t, Neutral *neutral);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Get the number of neutral species.
 *
 * @param neutral The neutral data.
 *
 * @return The number of neutral species.
 */
size_t Neutral_get_n_species(Neutral *neutral);

GPU_DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Get charge number for each neutral species.
 *
 * @param neutral The neutral data.
 *
 * @return The charge numbers.
 */
const int *Neutral_get_species_znum(Neutral *neutral);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Get atomic mass number for each neutral species.
 *
 * @param neutral The neutral data.
 *
 * @return The atomic mass numbers.
 */
const int *Neutral_get_species_anum(Neutral *neutral);
DECLARE_TARGET_END

#endif
