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
 * @param nspecies Number of neutral species.
 * @param nrho Number of rho grid points in the data.
 * @param anum Atomic mass number of the neutral species.
 * @param znum Charge number of the neutral species.
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
    NeutralRadial *neutral, size_t nspecies, size_t nrho, int anum[nspecies],
    int znum[nspecies], real rholim[2], real density[nspecies * nrho],
    real temperature[nspecies * nrho]);

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

GPU_DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate density of a neutral species.
 *
 * @param density Evaluated density [m^-3].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param i_species Index of the requested species.
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t NeutralRadial_eval_density(
    real density[1], real rho, size_t i_species, NeutralRadial *neutral);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate temperature of a neutral species.
 *
 * @param temperature Evaluated temperature [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param i_species Index of the requested species.
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t NeutralRadial_eval_temperature(
    real temperature[1], real rho, size_t i_species, NeutralRadial *neutral);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate density and temperature for all neutral species.
 *
 * @param density Evaluated density [m^-3].
 * @param temperature Evaluated temperature [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t NeutralRadial_eval_nT(
    real *density, real *temperature, real rho, NeutralRadial *neutral);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

#endif
