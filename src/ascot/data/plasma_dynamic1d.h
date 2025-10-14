/**
 * @file plasma_dynamic1d.h
 * @brief Header file for plasma_1Dt.c
 */
#ifndef PLASMA_DYNAMIC1D_H
#define PLASMA_DYNAMIC1D_H

#include "defines.h"
#include "parallel.h"
#include "plasma.h"

/**
 * Initialize the linearly interpolated radial time-dependent plasma data.
 *
 * @param plasma The struct to initialize.
 * @param nrho Number of rho grid points in the data.
 * @param ntime Number of time grid points in the data
 * @param nion Number of ion species.
 * @param rho Grid in rho in which data is tabulated [1].
 *        No need to be uniform.
 * @param time Grid in time in which data is tabulated [s].
 *        No need to be uniform.
 * @param anum Atomic mass number of the ion species.
 * @param znum Charge number of the ion species.
 * @param mass Mass of the ion species [kg].
 * @param charge Charge of the ion species [C].
 * @param Te Electron temperature [J].
 *        Layout is (rhoi, tj) = [j*nrho + i] (C order).
 * @param Ti Temperature of the ion species [J].
 *        Layout is (rhoi, tj) = [j*nrho + i] (C order).
 * @param ne Electron density [m^-3].
 *        Layout is (rhoi, tj) = [j*nrho + i] (C order).
 * @param ni Density of the ion species [m^-3].
 *        Layout is (ion, rhoi, tj) = [ion*nrho*ntime + j*nrho + i] (C order).
 * @param vtor Toroidal rotation of the whole plasma [rad/s].
 *        Layout is (rhoi, tj) = [j*nrho + i] (C order).
 *
 * @return Zero if the initialization succeeded.
 */
int PlasmaDynamic1D_init(
    PlasmaDynamic1D *plasma, size_t nrho, size_t ntime, size_t nion,
    int anum[nion], int znum[nion], real mass[nion], real charge[nion],
    real rho[nrho], real time[ntime], real Te[nrho * ntime],
    real Ti[nrho * ntime], real ne[nrho * ntime], real ni[nrho * ntime * nion],
    real vtor[nrho * ntime]);

/**
 * Free allocated resources.
 *
 * @param plasma The struct whose fields are deallocated.
 */
void PlasmaDynamic1D_free(PlasmaDynamic1D *plasma);

/**
 * Offload data to the accelerator.
 *
 * @param plasma The struct to offload.
 */
void PlasmaDynamic1D_offload(PlasmaDynamic1D *plasma);

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate temperature of a plasma species.
 *
 * @param temperature Evaluated temperature [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param t Time coordinate of the query point [s].
 * @param i_species Index of the requested species.
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaDynamic1D_eval_temperature(
    real temperature[1], real rho, real t, size_t i_species,
    PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate density of a plasma species.
 *
 * @param density Evaluated density [m^-3].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param t Time coordinate of the query point [s].
 * @param i_species Index of the requested species.
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaDynamic1D_eval_density(
    real density[1], real rho, real t, size_t i_species,
    PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma density and temperature for all species.
 *
 * @param density Evaluated density (electrons first followed by ions) [m^-3].
 * @param temperature Evaluated temperature (electrons first followed by ions)
 *        [J].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param t Time coordinate of the query point [s].
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaDynamic1D_eval_nT(
    real *density, real *temperature, real rho, real t,
    PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma flow along the field lines (same for all species).
 *
 * @param vflow Evaluated flow value [m/s].
 * @param rho Normalized poloidal flux coordinate of the query point [1].
 * @param r R coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaDynamic1D_eval_flow(
    real vflow[1], real rho, real t, real r, PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

#endif
