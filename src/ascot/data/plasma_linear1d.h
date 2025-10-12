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
 * Initialize 1D plasma data and check inputs
 *
 * @param data pointer to the data struct
 * @param nrho number of rho grid points in the data
 * @param nion number of ion species
 * @param rho rho grid in which data is tabulated [1]
 * @param anum mass number of ions present in the plasma
 * @param znum charge number of ions present in the plasma
 * @param mass mass of ions present in the plasma [kg]
 * @param charge charge of ions present in the plasma [C]
 * @param Te electron temperature [J]
 * @param Ti ion temperature [J]
 * @param ne electron density [m^-3]
 * @param ni density of ion species [m^-3]
 * @param vtor toroidal rotation [rad/s]
 *
 * @return zero if initialization succeeded
 */
int PlasmaLinear1D_init(
    PlasmaLinear1D *plasma, size_t nrho, size_t nion, int anum[nion], int znum[nion],
    real mass[nion], real charge[nion], real rho[nrho], real Te[nrho],
    real Ti[nrho], real ne[nrho], real ni[nrho * nion], real vtor[nrho]);

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
 * Evaluate plasma temperature
 *
 * This function evaluates the temperature of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param temp pointer to where evaluated temperature is stored [J]
 * @param rho radial coordinate [1]
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
err_t PlasmaLinear1D_eval_temp(
    real *dens, real rho, size_t species, PlasmaLinear1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma density
 *
 * This function evaluates the density of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param dens pointer to where evaluated density is stored [m^-3]
 * @param rho radial coordinate [1]
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
err_t PlasmaLinear1D_eval_dens(
    real *temp, real rho, size_t species, PlasmaLinear1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evaluate plasma density and temperature for all species
 *
 * This function evaluates the density and temperature of all plasma species at
 * the given radial coordinate using linear interpolation.
 *
 * @param dens pointer to where interpolated densities are stored [m^-3]
 * @param temp pointer to where interpolated temperatures are stored [J]
 * @param rho radial coordinate [1]
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
err_t PlasmaLinear1D_eval_densandtemp(
    real *dens, real *temp, real rho, PlasmaLinear1D *plasma);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
/**
 * Evalate plasma flow along the field lines
 *
 * @param vflow pointer where the flow value is stored [m/s]
 * @param rho particle rho coordinate [1]
 * @param r particle R coordinate [m]
 * @param pls_data pointer to plasma data
 */
err_t PlasmaLinear1D_eval_flow(
    real vflow[1], real rho, real r, PlasmaLinear1D *plasma);
DECLARE_TARGET_END

#endif
