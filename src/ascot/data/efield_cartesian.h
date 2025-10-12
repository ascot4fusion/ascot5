/**
 * @file efield_cartesian.h
 * Cartesian electric field implementation.
 *
 * This electric field is constant everywhere but the field vector is defined
 * on a Cartesian basis (so the vector rotates in the cylindrical basis). This
 * field is mainly used for testing the orbit integration.
 */
#ifndef EFIELD_CARTESIAN_H
#define EFIELD_CARTESIAN_H

#include "bfield.h"
#include "defines.h"
#include "efield.h"
#include "parallel.h"

/**
 * Initialize the Cartesian electric field.
 *
 * Nothing is allocated; only the values of ``exyz`` are copied.
 *
 * @param efield The struct to initialize.
 * @param exyz Electric field vector on cartesian basis [V/m].
 *
 * @return Zero if the initialization succeeded.
 */
int EfieldCartesian_init(EfieldCartesian *efield, real exyz[3]);

/**
 * Free allocated resources.
 *
 * Does nothing as no resources were dynamically allocated.
 *
 * @param efield The struct whose fields are deallocated.
 */
void EfieldCartesian_free(EfieldCartesian *efield);

/**
 * Offload data to the accelerator.
 *
 * @param efield The struct to offload.
 */
void EfieldCartesian_offload(EfieldCartesian *efield);

GPU_DECLARE_TARGET_SIMD_UNIFORM(efield)
/**
 * Evaluate electric field vector.
 *
 * The electric field is constant everywhere, so the stored field vector is
 * converted to cylindrical coordinates and returned.
 *
 * @param e Evaluated electric field [V/m].
 *        Layout: [er, ephi, ez].
 * @param phi phi coordinate of the query point [rad].
 * @param efield The electric field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t EfieldCartesian_eval_e(real e[3], real phi, EfieldCartesian *efield);
DECLARE_TARGET_END

#endif
