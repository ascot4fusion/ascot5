/**
 * @file efield_potential1d.h
 * Radial electric field evaluated from the gradient of a 1D potential.
 */
#ifndef EFIELD_POTENTIAL1D_H
#define EFIELD_POTENTIAL1D_H
#include "bfield.h"
#include "defines.h"
#include "efield.h"
#include "parallel.h"
#include "utils/interp.h"

/**
 * Initialize the 1D potential electric field.
 *
 * Allocates the linear interpolant used to evaluate the electric field.
 *
 * @param efield The struct to initialize.
 * @param nrho Number of points in the rho grid.
 * @param rholim Range of the uniform rho grid [1].
 * @param dvdrho The gradient of the electric field potential with respect to
 *        rho [V].
 *
 * @return Zero if the initialization succeeded.
 */
int EfieldPotential1D_init(
    EfieldPotential1D *efield, size_t nrho, real rholim[2], real dvdrho[nrho]);

/**
 * Free allocated resources.
 *
 * @param efield The struct whose fields are deallocated.
 */
void EfieldPotential1D_free(EfieldPotential1D *efield);

/**
 * Offload data to the accelerator.
 *
 * @param efield The struct to offload.
 */
void EfieldPotential1D_offload(EfieldPotential1D *efield);

GPU_DECLARE_TARGET_SIMD_UNIFORM(efield, bfield)
/**
 * Evaluate electric field vector.
 *
 * This function first evaluates the gradient of rho, which is needed to convert
 * the gradient of the electric potential (with respect to rho) to V/m.
 *
 * @param e Evaluated electric field [V/m].
 *        Layout: [er, ephi, ez].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param efield The electric field data.
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t EfieldPotential1D_eval_e(
    real e[3], real r, real phi, real z, real t, EfieldPotential1D *efield,
    Bfield *bfield);
DECLARE_TARGET_END

#endif
