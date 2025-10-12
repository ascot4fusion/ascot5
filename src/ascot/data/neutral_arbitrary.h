/**
 * @file neutral_arbitrary.h
 * 3D neutral data with trilinear interpolation.
 */
#ifndef NEUTRAL_ARBITRARY_H
#define NEUTRAL_ARBITRARY_H

#include "defines.h"
#include "neutral.h"
#include "parallel.h"

/**
 * Initialize neutral data.
 *
 * @param data pointer to data struct.
 * @param n_r number of r grid points in the data.
 * @param r_min minimum r coordinate in the grid in the data [m].
 * @param r_max maximum r coordinate in the grid in the data [m].
 * @param n_phi number of phi grid points in the data.
 * @param phi_min minimum phi coordinate in the grid in the data [rad].
 * @param phi_max maximum phi coordinate in the grid in the data [rad].
 * @param n_z number of z grid points in the data.
 * @param z_min minimum z coordinate in the grid in the data [m].
 * @param z_max maximum z coordinate in the grid in the data [m].
 * @param n_species number of neutral species.
 * @param anum neutral species mass number.
 * @param znum neutral species charge number.
 * @param density neutral species-wise density [m^-3].
 * @param temperature neutral species-wise temperature [J].
 */
int NeutralArbitrary_init(
    NeutralArbitrary *neutral, size_t n, size_t nr, size_t nphi, size_t nz,
    real rlim[2], real philim[2], real zlim[2], int *anum, int *znum,
    real *density, real *temperature);

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void NeutralArbitrary_free(NeutralArbitrary *neutral);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void NeutralArbitrary_offload(NeutralArbitrary *neutral);

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density at the given coordinates using
 * trilinear interpolation on the 3D neutral density data.
 *
 * @param n0 pointer where neutral density is stored [m^-3]
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param ndata pointer to neutral density data struct
 *
 * @return zero if evaluation succeeded
 */
DECLARE_TARGET_SIMD_UNIFORM(neutral)
err_t NeutralArbitrary_eval_n0(real *n0, real r, real phi, real z, NeutralArbitrary *neutral);

/**
 * @brief Evaluate neutral temperature
 *
 * This function evaluates the neutral temperature at the given coordinates
 * using trilinear interpolation on the 3D neutral temperature data.
 *
 * @param t0 pointer where neutral temperature is stored [J]
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 *
 * @return zero if evaluation succeeded
 */
DECLARE_TARGET_SIMD_UNIFORM(neutral)
err_t NeutralArbitrary_eval_T0(real *T0, real r, real phi, real z, NeutralArbitrary *neutral);

#endif
