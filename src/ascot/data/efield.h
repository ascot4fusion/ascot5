/**
 * @file efield.h
 * Electric field interface.
 *
 * Provides functions and datatypes for electric field evaluation.
 */
#ifndef EFIELD_H
#define EFIELD_H

#include "bfield.h"
#include "defines.h"
#include "efield_data.h"
#include "parallel.h"

/**
 * Free allocated resources.
 *
 * @param efield The struct whose fields are deallocated.
 */
void Efield_free(Efield *efield);

/**
 * Offload data to the accelerator.
 *
 * @param efield The struct to offload.
 */
void Efield_offload(Efield *efield);

GPU_DECLARE_TARGET_SIMD_UNIFORM(efield, bfield)
/**
 * Evaluate electric field vector.
 *
 * If electric field input is not specified, instead of raising error a vector
 * that is zero is returned.
 *
 * @param e Evaluated electric field [V/m].
 *        Layout: [er, ephi, ez].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param t Time coordinate of the query point [s].
 * @param efield The electric field data.
 * @param bfield The magnetic field data.
 *        This is required as the electric field may be a function of normalized
 *        poloidal flux.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t Efield_eval_e(
    real e[3], real r, real phi, real z, real t, Efield *efield,
    Bfield *bfield);
DECLARE_TARGET_END

#endif
