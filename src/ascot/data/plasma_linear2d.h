/**
 * @file plasma_linear2D.h
 * Header file for PlasmaLinear2D.c
 */
#ifndef PlasmaLinear2D_H
#define PlasmaLinear2D_H

#include "defines.h"
#include "parallel.h"
#include "plasma.h"

/**
 * Initialize the linearly interpolated axisymmetric plasma data.
 *
 * @param plasma The struct to initialize.
 * @param nr Number of radial grid points in the data.
 * @param nz Number of axial grid points in the data.
 * @param nion Number of ion species.
 * @param rmin Grid in rho in which data is tabulated [1].
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
int PlasmaLinear2D_init(
    PlasmaLinear2D *plasma, size_t nr, size_t nz, size_t nion, real rlim[2],
    real zlim[2], int anum[nion], int znum[nion], real mass[nion],
    real charge[nion], real Te[nr * nz], real Ti[nr * nz], real ne[nr * nz],
    real ni[nr * nz * nion], real vtor[nr * nz]);

/**
 * Free allocated resources.
 *
 * @param plasma The struct whose fields are deallocated.
 */
void PlasmaLinear2D_free(PlasmaLinear2D *plasma);

/**
 * Offload data to the accelerator.
 *
 * @param plasma The struct to offload.
 */
void PlasmaLinear2D_offload(PlasmaLinear2D *plasma);

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate temperature of a plasma species.
 *
 * @param temperature Evaluated temperature [J].
 * @param r Radial coordinate of the query point [m].
 * @param z Axial coordinate of the query point [m].
 * @param i_species Index of the requested species.
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear2D_eval_temperature(
    real temperature[1], real r, real z, size_t i_species,
    PlasmaLinear2D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate density of a plasma species.
 *
 * @param density Evaluated density [m^-3].
 * @param r Radial coordinate of the query point [m].
 * @param z Axial coordinate of the query point [m].
 * @param i_species Index of the requested species.
 *        Zero for electrons, then 1 for the first ion, etc. in the same order
 *        as they are listed in the plasma data.
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear2D_eval_density(
    real density[1], real r, real z, size_t i_species, PlasmaLinear2D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma density and temperature for all species.
 *
 * @param density Evaluated density (electrons first followed by ions) [m^-3].
 * @param temperature Evaluated temperature (electrons first followed by ions)
 *        [J].
 * @param r Radial coordinate of the query point [m].
 * @param z Axial coordinate of the query point [m].
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear2D_eval_nT(
    real *density, real *temperature, real r, real z, PlasmaLinear2D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma flow along the field lines (same for all species).
 *
 * @param vflow Evaluated flow value [m/s].
 * @param r Radial coordinate of the query point [m].
 * @param z Axial coordinate of the query point [m].
 * @param plasma The plasma data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t PlasmaLinear2D_eval_flow(
    real vflow[1], real r, real z, PlasmaLinear2D *plasma);
DECLARE_TARGET_END

#endif