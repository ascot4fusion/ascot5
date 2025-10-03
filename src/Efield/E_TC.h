/**
 * @brief Header file for E_TC.c
 *
 * Contains declaration of E_TC_field_offload_data and E_TC_field_data structs.
 */
#ifndef E_TC_H
#define E_TC_H

#include "offload.h"
#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "E_field.h"


/**
 * @brief Initialize electric field data and check inputs
 *
 * @param data pointer to the data struct
 * @param exyz electric field vector on cartesian basis [V/m]
 *
 * @return Zero to indicate initialization succeeded
 */
int EfieldCartesian_init(EfieldCartesian* efield, real exyz[3]);


/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void EfieldCartesian_free(EfieldCartesian* efield);


/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void EfieldCartesian_offload(EfieldCartesian* efield);


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
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Zero to indicate success
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(efield, bfield)
a5err EfieldCartesian_eval_e(
    real e[3], real r, real phi, real z, EfieldCartesian* efield,
    B_field_data* bfield
);
DECLARE_TARGET_END
#endif
