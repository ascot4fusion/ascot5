/**
 * @file N0_1D.h
 * 1D neutral data with linear interpolation.
 */
#ifndef N0_1D_H
#define N0_1D_H
#include "defines.h"
#include "parallel.h"
#include "neutral.h"

/**
 * Initialize data.
 *
 * @param data pointer to data struct.
 * @param n_rho number of rho grid points in the data.
 * @param rho_min minimum rho coordinate in the grid in the data [1].
 * @param rho_max maximum rho coordinate in the grid in the data [1].
 * @param n_species number of neutral species.
 * @param anum neutral species mass number.
 * @param znum neutral species charge number.
 * @param density neutral species-wise density [m^-3].
 * @param temperature neutral species-wise temperature [J].
 *
 * @return Zero if initialization succeeded.
 */
int N0_1D_init(N0_1D_data* data, int n_rho, real rho_min, real rho_max,
               int n_species, int* anum, int* znum,
               real* density, real* temperature);

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void N0_1D_free(N0_1D_data* data);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void N0_1D_offload(N0_1D_data* data);

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density at the given coordinates using
 * linear interpolation on the 1D neutral density data.
 *
 * @param n0 pointer where neutral density is stored [m^-3]
 * @param rho normalized poloidal flux coordinate
 * @param ndata pointer to neutral data struct
 *
 * @return zero if evaluation succeeded
 */
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_1D_eval_n0(real* n0, real rho, N0_1D_data* ndata);

/**
 * @brief Evaluate neutral temperature
 *
 * This function evaluates the neutral temperature at the given coordinates
 * using linear interpolation on the 1D neutral temperature data.
 *
 * @param t0 pointer where neutral temperature is stored [J]
 * @param rho normalized poloidal flux coordinate
 *
 * @return zero if evaluation succeeded
 */
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_1D_eval_t0(real* t0, real rho, N0_1D_data* ndata);

/**
 * @brief Return number of neutral species
 *
 * @param ndata pointer to neutral data struct
 *
 * @return number of neutral species
 */
DECLARE_TARGET_SIMD_UNIFORM(ndata)
int N0_1D_get_n_species(N0_1D_data* ndata);

#endif
