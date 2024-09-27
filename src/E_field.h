/**
 * @file E_field.h
 * @brief Header file for E_field.c
 *
 * Contains a list declaring all E_field_types, and declaration of
 * E_field_offload_data and E_field_data structs.
 */
#ifndef E_FIELD_H
#define E_FIELD_H

#include "offload_acc_omp.h"
#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "Efield/E_TC.h"
#include "Efield/E_1DS.h"

/**
 * @brief Electric field types
 *
 * Electric field types are used in the electric field interface (E_field.c) to
 * direct function calls to correct electric field instances. Each electric
 * field instance must have a corresponding type.
 */
typedef enum E_field_type {
    E_field_type_TC, /**< Trivial Cartesian electric field */
    E_field_type_1DS /**< Spline-interpolated radial electric field */
} E_field_type;

/**
 * @brief Electric field simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct {
    E_field_type type; /**< Electric field type wrapped by this struct */
    E_TC_data ETC;     /**< TC field or NULL if not active             */
    E_1DS_data E1DS;   /**< 1DS field or NULL if not active            */
} E_field_data;

void E_field_free(E_field_data* Edata);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Edata, Bdata)
a5err E_field_eval_E(real E[3], real r, real phi, real z, real t,
                     E_field_data* Edata, B_field_data* Bdata);
DECLARE_TARGET_END

#endif
