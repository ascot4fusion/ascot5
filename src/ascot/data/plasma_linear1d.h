/**
 * @file plasma_linear1D.h
 * Header file for PlasmaLinear1D.c
 */
#ifndef PlasmaLinear1D_H
#define PlasmaLinear1D_H

#include "defines.h"
#include "parallel.h"
#include "plasma.h"

/**
 * Initialize the linearly interpolated radial plasma data.
 *
 * @param plasma The struct to initialize.
 * @param nrho Number of rho grid points in the data.
 * @param nion Number of ion species.
 * @param rho Grid in rho in which data is tabulated [1].
 *        No need to be uniform.
 * @param anum Atomic mass number of the ion species.
 * @param znum Charge number of the ion species.
 * @param mass Mass of the ion species [kg].
 * @param charge Charge of the ion species [C].
 * @param Te Electron temperature [J].
 * @param Ti Temperature of the ion species [J].
 * @param ne Electron density [m^-3].
 * @param ni Density of the ion species [m^-3].
 *        Layout is (ion, rhoi) = [ion*nrho + i] (C order).
 * @param vtor Toroidal rotation of the whole plasma [rad/s].
 *
 * @return Zero if the initialization succeeded.
 */
int PlasmaLinear1D_init(
    PlasmaLinear1D *plasma, size_t nrho, size_t nion, int anum[nion],
    int znum[nion], real mass[nion], real charge[nion], real rho[nrho],
    real Te[nrho], real Ti[nrho], real ne[nrho], real ni[nrho * nion],
    real vtor[nrho]);

/**
 * Free allocated resources.
 *
 * @param plasma The struct whose fields are deallocated.
 */
void PlasmaLinear1D_free(PlasmaLinear1D *plasma);

/**
 * Offload data to the accelerator.
 *
 * @param plasma The struct to offload.
 */
void PlasmaLinear1D_offload(PlasmaLinear1D *plasma);

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate temperature of a plasma species.
 *
 * @param temperature Evaluated temperature [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param i_species Index of the requested species.
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear1D_eval_temperature(
    real temperature[1], real rho, size_t i_species, PlasmaLinear1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate density of a plasma species.
 *
 * @param density Evaluated density [m^-3].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param i_species Index of the requested species.
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear1D_eval_density(
    real density[1], real rho, size_t i_species, PlasmaLinear1D *plasma);
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
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear1D_eval_nT(
    real *density, real *temperature, real rho, PlasmaLinear1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma flow along the field lines (same for all species).
 *
 * @param vflow Evaluated flow value [m/s].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear1D_eval_flow(
    real vflow[1], real rho, real r, PlasmaLinear1D *plasma);
DECLARE_TARGET_END

#endif
