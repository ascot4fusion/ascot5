/**
 * @file neutral_arbitrary.h
 * Arbitrary (3D) neutral data interpolated linearly.
 */
#ifndef NEUTRAL_ARBITRARY_H
#define NEUTRAL_ARBITRARY_H

#include "defines.h"
#include "neutral.h"
#include "parallel.h"

/**
 * Initialize the linearly interpolated cylindrical neutral data.
 *
 * @param neutral The struct to initialize.
 * @param n Number of neutral species.
 * @param nr Number of R grid points in the data.
 * @param nphi Number of phi grid points in the data.
 * @param nz Number of z grid points in the data.
 * @param rlim Limits of the uniform R grid [m].
 * @param philim Limits of the uniform phi grid [rad].
 * @param zlim Limits of the uniform z grid [m].
 * @param density Density of the neutral species (in same order as they are
 *        listed in the plasma input) [m^-3].
 *        Layout is (neutral, rhoi) = [neutral*nrho + i] (C order).
 * @param temperature Temperature of the neutral species (in same order as they
 *        are listed in the plasma input) [J].
 *        Layout is (neutral, rhoi) = [neutral*nrho + i] (C order).
 *
 * @return Zero if initialization succeeded.
 */
int NeutralArbitrary_init(
    NeutralArbitrary *neutral, size_t n, size_t nr, size_t nphi, size_t nz,
    real rlim[2], real philim[2], real zlim[2],
    real density[n * nr * nphi * nz], real temperature[n * nr * nphi * nz]);

/**
 * Free allocated resources.
 *
 * @param neutral The struct whose fields are deallocated.
 */
void NeutralArbitrary_free(NeutralArbitrary *neutral);

/**
 * Offload data to the accelerator.
 *
 * @param neutral The struct to offload.
 */
void NeutralArbitrary_offload(NeutralArbitrary *neutral);

DECLARE_TARGET_SIMD_UNIFORM(neutral)
/**
 * Evaluate density of each neutral species.
 *
 * The evaluated values are in the same order as ions in the plasma data.
 *
 * @param density Evaluated density [m^-3].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t NeutralArbitrary_eval_density(
    real *density, real r, real phi, real z, NeutralArbitrary *neutral);

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
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param neutral The neutral data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t NeutralArbitrary_eval_temperature(
    real *temperature, real r, real phi, real z, NeutralArbitrary *neutral);

#endif
