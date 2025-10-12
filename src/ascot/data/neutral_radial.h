/**
 * @file neutral_radial.h
 * 1D neutral data with linear interpolation.
 */
#ifndef NEUTRAL_RADIAL_H
#define NEUTRAL_RADIAL_H

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
int NeutralRadial_init(
    NeutralRadial *neutral, size_t n, size_t nrho, int *anum, int *znum,
    real rholim[2], real *density, real *temperature);

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void NeutralRadial_free(NeutralRadial* neutral);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void NeutralRadial_offload(NeutralRadial* neutral);

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
DECLARE_TARGET_SIMD_UNIFORM(neutral)
err_t NeutralRadial_eval_n0(real* n0, real rho, NeutralRadial* neutral);

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
DECLARE_TARGET_SIMD_UNIFORM(neutral)
err_t NeutralRadial_eval_T0(real* T0, real rho, NeutralRadial* neutral);

#endif
