/**
 * @file plasma_dynamic1d.h
 * @brief Header file for plasma_1Dt.c
 */
#ifndef PLASMA_DYNAMIC1D_H
#define PLASMA_DYNAMIC1D_H
#include "defines.h"
#include "errors.h"
#include "parallel.h"
#include "plasma.h"

/**
 * @brief Initialize 1Dt plasma data and check inputs
 *
 * @param data pointer to the data struct
 * @param nrho number of rho grid points in the data
 * @param nrho number of time grid points in the data
 * @param nion number of ion species
 * @param rho rho grid points in which data is tabulated [1]
 * @param time time grid points in which data is tabulated [s]
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
int PlasmaDynamic1D_init(
    PlasmaDynamic1D *plasma, int nrho, int ntime, int nion, int *anum,
    int *znum, real *rho, real *time, real *mass, real *charge, real *Te,
    real *Ti, real *ne, real *ni, real *vtor);

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void PlasmaDynamic1D_free(PlasmaDynamic1D *plasma);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void PlasmaDynamic1D_offload(PlasmaDynamic1D *plasma);

/**
 * @brief Evaluate plasma temperature
 *
 * This function evaluates the temperature of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param temp pointer to where evaluated temperature is stored [J]
 * @param rho radial coordinate [1]
 * @param t time instant [s]
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err PlasmaDynamic1D_eval_temp(
    real *dens, real rho, real t, int species, PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

/**
 * @brief Evaluate plasma density
 *
 * This function evaluates the density of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param dens pointer to where evaluated density is stored [m^-3]
 * @param rho radial coordinate [1]
 * @param t time instant [s]
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err PlasmaDynamic1D_eval_dens(
    real *temp, real rho, real t, int species, PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 * This function evaluates the density and temperature of all plasma species at
 * the given radial coordinate using linear interpolation.
 *
 * @param dens pointer to where interpolated densities are stored [m^-3]
 * @param temp pointer to where interpolated temperatures are stored [J]
 * @param rho radial coordinate [1]
 * @param t time instant [s]
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err PlasmaDynamic1D_eval_densandtemp(
    real *dens, real *temp, real rho, real t, PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

/**
 * @brief Evalate plasma flow along the field lines
 *
 * @param vflow pointer where the flow value is stored [m/s]
 * @param rho particle rho coordinate [1]
 * @param t particle time coordinate [s]
 * @param r particle R coordinate [m]
 * @param pls_data pointer to plasma data
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err PlasmaDynamic1D_eval_flow(
    real *vflow, real rho, real t, real r, PlasmaDynamic1D *plasma);
DECLARE_TARGET_END

#endif
