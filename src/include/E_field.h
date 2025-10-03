/**
 * Contains a list declaring all E_field_types, and declaration of
 * E_field_offload_data and E_field_data structs.
 */
#ifndef E_FIELD_H
#define E_FIELD_H

#include "B_field.h"
#include "ascot5.h"
#include "error.h"
#include "offload.h"

/**
 * Electric field types
 *
 * Electric field types are used in the electric field interface (E_field.c) to
 * direct function calls to correct electric field instances. Each electric
 * field instance must have a corresponding type.
 */
typedef enum E_field_type
{
    E_field_type_cartesian,  /**< Trivial Cartesian electric field */
    E_field_type_potential1d /**< Spline-interpolated radial electric field */
} E_field_type;

/**
 * @brief Trivial Cartesian electric field simulation data
 */
typedef struct {
    real exyz[3]; /**< Pointer to array holding constant [E_x, E_y, E_z]
                       values [V/m]                                           */
} EfieldCartesian;

/**
 * @brief 1D spline electric field parameters on the target
 */
typedef struct {
    real reff;         /**< Effective minor radius [m] */
    interp1D_data dv;  /**< dV_drho 1D linear interpolation struct */
} EfieldPotential1D;

/**
 * @brief Electric field simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct
{
    EfieldCartesian *cartesian;     /**< TC field or NULL if not active     */
    EfieldPotential1D *potential1d; /**< 1DS field or NULL if not active */
    E_field_type type; /**< Electric field type wrapped by this struct */
} E_field_data;

void E_field_free(E_field_data *data);

void E_field_offload(E_field_data *data);

GPU_DECLARE_TARGET_SIMD_UNIFORM(efield, bfield)
a5err E_field_eval_E(
    real e[3], real r, real phi, real z, real t, E_field_data *efield,
    B_field_data *bfield);
DECLARE_TARGET_END

#endif
