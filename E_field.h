/**
 * @file E_field.h
 * @brief Header file for E_field.c
 *
 * Contains a list declaring all E_field_types, and declaration of
 * E_field_offload_data and E_field_data structs.
 */
#ifndef E_FIELD_H
#define E_FIELD_H

#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "Efield/E_TC.h"
#include "Efield/E_1DS.h"
#include "Efield/E_3D.h"
#include "Efield/E_3DS.h"
#include "Efield/E_3DST.h"

/**
 * @brief Electric field types
 *
 * Electric field types are used in the electric field interface (E_field.c) to
 * direct function calls to correct electric field instances. Each electric
 * field instance must have a corresponding type.
 */
typedef enum E_field_type {
    E_field_type_TC,   /**< Trivial Cartesian electric field          */
    E_field_type_1DS,  /**< Spline-interpolated radial electric field */
    E_field_type_3D,   /**< Linear-interpolated 3D electric field     */
    E_field_type_3DS,  /**< Spline-interpolated 3D electric field     */
    E_field_type_3DST  /**< 3D time-dependent electric field          */
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
    E_3D_offload_data E3D;    /**< 3D fiel or NULL if not active              */
    E_3DS_offload_data E3DS;  /**< 3DS fiel or NULL if not active             */
    E_3DST_offload_data E3DST;/**< 3DST fiel or NULL if not active            */
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
    E_field_type type; /**< Electric field type wrapped by this struct  */
    E_TC_data ETC;     /**< TC field or NULL if not active              */
    E_1DS_data E1DS;   /**< 1DS field or NULL if not active             */
    E_3D_data E3D;     /**< 3D field or NULL if not active              */
    E_3DS_data E3DS;   /**< 3DS field or NULL if not active             */
    E_3DST_data E3DST;  /**< 3DST field or NULL if not active             */

} E_field_data;

int E_field_init_offload(E_field_offload_data* offload_data,
                         real** offload_array);
void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
int E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
                 real* offload_array);
#pragma omp declare simd uniform(Edata, Bdata)
a5err E_field_eval_E(real* E, real r, real phi, real z, real t,
                     E_field_data* Edata, B_field_data* Bdata);
#pragma omp end declare target

#endif
