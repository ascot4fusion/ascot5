/**
 * @file N0_3D.h
 * 3D neutral data with trilinear interpolation.
 */
#ifndef N0_3D_H
#define N0_3D_H
#include "ascot5.h"
#include "error.h"
#include "neutral.h"
#include "offload.h"

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
int N0_3D_init(
    N0_3D_data *data, int n_r, real r_min, real r_max, int n_phi, real phi_min,
    real phi_max, int n_z, real z_min, real z_max, int n_species, int *anum,
    int *znum, real *density, real *temperature);

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void N0_3D_free(N0_3D_data *data);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void N0_3D_offload(N0_3D_data *data);

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
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_3D_eval_n0(real *n0, real r, real phi, real z, N0_3D_data *ndata);

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
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_3D_eval_t0(real *t0, real r, real phi, real z, N0_3D_data *ndata);

/**
 * @brief Return number of neutral species
 *
 * @param ndata pointer to neutral data struct
 *
 * @return number of neutral species
 */
DECLARE_TARGET_SIMD_UNIFORM(ndata)
int N0_3D_get_n_species(N0_3D_data *ndata);
#endif
