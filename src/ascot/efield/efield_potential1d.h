/**
 * @file efield_potential1d.h
 * Contains declaration of E_1DS_field_offload_data and E_1DS_field_data
 * structs.
 */
#ifndef EFIELD_POTENTIAL1D_H
#define EFIELD_POTENTIAL1D_H
#include "defines.h"
#include "bfield.h"
#include "efield.h"
#include "errors.h"
#include "interp.h"
#include "parallel.h"

/**
 * @brief Initialize 1DS electric field data
 *
 * @param data pointer to the data struct
 * @param nrho number of points in the rho grid
 * @param rhomin minimum rho value in the grid
 * @param rhomax maximum rho value in the grid
 * @param reff effective minor radius
 * @param dvdrho gradient of the potential
 *
 * @return zero if initialization succeeded
 */
int EfieldPotential1D_init(
    EfieldPotential1D *efield, int nrho, real reff, real rholim[2],
    real dvdrho[nrho]);

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void EfieldPotential1D_free(EfieldPotential1D *efield);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void EfieldPotential1D_offload(EfieldPotential1D *efield);

/**
 * @brief Evaluate 1D spline radial electric field
 *
 * This function evaluates the 1D spline potential gradient of the plasma at the
 * given radial coordinate using linear interpolation, and then calculates the
 * radial electric field by multiplying that with the rho-gradient. Gradient of
 * rho is obtained via magnetic field module.
 *
 * @param E array where the electric field will be stored (E_r -> E[1],
 *        E_phi -> E[1], E_z -> E[2])
 * @param r R-coordiante [m]
 * @param phi phi-coordinate [rad]
 * @param z z-coordiante [m]
 * @param Edata pointer to electric field data
 * @param Bdata pointer to magnetic field data
 *
 * @return zero if evaluation succeeded
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(efield, bfield)
a5err EfieldPotential1D_eval_e(
    real e[3], real r, real phi, real z, EfieldPotential1D *efield,
    Bfield *bfield);
DECLARE_TARGET_END

#endif
