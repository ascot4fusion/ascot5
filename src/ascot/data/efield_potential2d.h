/**
 * @file efield_potential2d.h
 * Electric field evaluated from 2D potential.
 */
#ifndef EFIELD_POTENTIAL2D_H
#define EFIELD_POTENTIAL2D_H
#include "defines.h"
#include "efield.h"
#include "parallel.h"

/**
 * Initialize the 2D potential electric field.
 *
 * Allocates the spline interpolant used to evaluate the electric field.
 *
 * @param efield The struct to initialize.
 * @param nrho Number of points in the radial grid.
 * @param nrho Number of points in the axial grid.
 * @param rlim Range of the uniform radial grid [m].
 * @param zlim Range of the uniform axial grid [m].
 * @param vpot The electric field potential in the grid points.
 *        Layout: (Ri, zj) = [i*nz + j] (C order).
 *
 * @return Zero if the initialization succeeded.
 */
int EfieldPotential2D_init(
    EfieldPotential2D *efield, size_t nr, size_t nz, real rlim[2],
    real zlim[2], real vpot[nr*nz]);

/**
 * Free allocated resources.
 *
 * @param efield The struct whose fields are deallocated.
 */
void EfieldPotential2D_free(EfieldPotential2D *efield);

/**
 * Offload data to the accelerator.
 *
 * @param efield The struct to offload.
 */
void EfieldPotential2D_offload(EfieldPotential2D *efield);

GPU_DECLARE_TARGET_SIMD_UNIFORM(efield)
/**
 * Evaluate electric field vector.
 *
 * @param e Evaluated electric field [V/m].
 *        Layout: [er, ephi, ez].
 * @param r R coordinate of the query point [m].
 * @param z z coordinate of the query point [m].
 * @param efield The electric field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t EfieldPotential2D_eval_e(
    real e[3], real r, real z, EfieldPotential2D *efield);
DECLARE_TARGET_END

#endif
