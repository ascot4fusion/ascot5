/**
 * @file efield_cartesian.h
 *
 * Contains declaration of EFIELD_CARTESIAN_field_offload_data and EFIELD_CARTESIAN_field_data structs.
 */
#ifndef EFIELD_CARTESIAN_H
#define EFIELD_CARTESIAN_H

#include "defines.h"
#include "bfield.h"
#include "efield.h"
#include "errors.h"
#include "parallel.h"

/**
 * @brief Initialize electric field data and check inputs
 *
 * @param data pointer to the data struct
 * @param exyz electric field vector on cartesian basis [V/m]
 *
 * @return Zero to indicate initialization succeeded
 */
int EfieldCartesian_init(EfieldCartesian *efield, real exyz[3]);

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void EfieldCartesian_free(EfieldCartesian *efield);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void EfieldCartesian_offload(EfieldCartesian *efield);

/**
 * @brief Evaluate electric field
 *
 * Even though this module represents a Cartesian electric field, the returned
 * values are given in cylindrical coordinates.
 *
 * @param E pointer to array where electric field values are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param Edata pointer to magnetic field data struct
 *
 * @return Zero to indicate success
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(efield)
a5err EfieldCartesian_eval_e(real e[3], real phi, EfieldCartesian *efield);
DECLARE_TARGET_END

#endif
