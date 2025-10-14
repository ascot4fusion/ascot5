/**
 * @file plasma.h
 * Plasma interface.
 *
 * Provides functions and datatypes for getting plasma species information and
 * for interpolating background plasma density and temperature.
 */
#ifndef PLASMA_H
#define PLASMA_H

#include "defines.h"
#include "parallel.h"
#include "plasma_data.h"

/**
 * Free allocated resources.
 *
 * @param plasma The struct whose fields are deallocated.
 */
void Plasma_free(Plasma *plasma);

/**
 * Offload data to the accelerator.
 *
 * @param plasma The struct to offload.
 */
void Plasma_offload(Plasma *plasma);

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate temperature of a plasma species.
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
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Plasma_eval_temperature(
    real temperature[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate density of a plasma species.
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
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Plasma_eval_density(
    real density[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma density and temperature for all species.
 *
 * Depending on the implementation, either the normalized poloidal flux
 * coordinate or the cylindrical coordinate is used in the interpolation.
 *
 * @param density Evaluated density (electrons first followed by ions) [m^-3].
 * @param temperature Evaluated temperature (electrons first followed by ions)
 *        [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Plasma_eval_nT(
    real *density, real *temperature, real rho, real r, real phi, real z,
    real t, Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma flow along the field lines (same for all species).
 *
 * Depending on the implementation, either the normalized poloidal flux
 * coordinate or the cylindrical coordinate is used in the interpolation. In any
 * case the R coordinate is required to convert the rotation (rad/s) to flow
 * (m/s).
 *
 * @param vflow Evaluated flow value [m/s].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Plasma_eval_flow(
    real vflow[1], real rho, real r, real phi, real z, real t, Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Get the number of plasma species (including electrons).
 *
 * @param plasma The plasma data.
 *
 * @return The number of plasma species
 */
size_t Plasma_get_n_species(Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Get mass for all plasma species.
 *
 * @param plasma The plasma data.
 *
 * @return The mass (electrons first followed by ions) [kg].
 */
const real *Plasma_get_species_mass(Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Get charge for all plasma species.
 *
 * @param plasma The plasma data.
 *
 * @return The charge (electrons first followed by ions) [C].
 */
const real *Plasma_get_species_charge(Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Get charge number for each ion species.
 *
 * @param plasma The plasma data.
 *
 * @return The charge numbers.
 */
const int *Plasma_get_species_znum(Plasma *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Get atomic mass number for each ion species.
 *
 * @param plasma The plasma data.
 *
 * @return The atomic mass numbers.
 */
const int *Plasma_get_species_anum(Plasma *plasma);
DECLARE_TARGET_END

#endif
