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
 * @brief Electric field offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized in
 * E_field_init_offload().
 *
 * The intended usage is that only single offload data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    E_field_type type;        /**< Electric field type wrapped by this struct */
    E_TC_offload_data ETC;    /**< TC field or NULL if not active             */
    E_1DS_offload_data E1DS;  /**< 1DS field or NULL if not active            */
    int offload_array_length; /**< Allocated offload array length             */
} E_field_offload_data;

/**
 * @brief Electric field simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from the E_field_offload_data in E_field_init().
 *
 * The intended usage is that only single E_field_data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    E_field_type type; /**< Electric field type wrapped by this struct */
    E_TC_data ETC;     /**< TC field or NULL if not active             */
    E_1DS_data E1DS;   /**< 1DS field or NULL if not active            */
} E_field_data;

int E_field_init_offload(E_field_offload_data* offload_data,
                         real** offload_array);
void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array);

int E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
                 real* offload_array);
#ifndef GPU
DECLARE_TARGET_SIMD_UNIFORM(Edata, Bdata)
#else
DECLARE_TARGET
#endif
a5err E_field_eval_E(real E[3], real r, real phi, real z, real t,
                     E_field_data* Edata, B_field_data* Bdata);
DECLARE_TARGET_END

#endif
