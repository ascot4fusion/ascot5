/**
 * @file B_field.h
 * @brief Header file for B_field.c
 *
 * Contains a list declaring all B_field_types, and declaration of
 * B_field_offload_data and B_field_data structs.
 */
#ifndef B_FIELD_H
#define B_FIELD_H

#include "offload_acc_omp.h"
#include "ascot5.h"
#include "error.h"
#include "Bfield/B_GS.h"
#include "Bfield/B_2DS.h"
#include "Bfield/B_3DS.h"
#include "Bfield/B_STS.h"
#include "Bfield/B_TC.h"

/**
 * @brief Magnetic field types
 *
 * Magnetic field types are used in the magnetic field interface (B_field.c) to
 * direct function calls to correct magnetic field instances. Each magnetic
 * field instance must have a corresponding type.
 */
typedef enum B_field_type {
    B_field_type_GS,  /**< Analytic magnetic field                          */
    B_field_type_2DS, /**< Spline-interpolated axisymmetric  magnetic field */
    B_field_type_3DS, /**< Spline-interpolated 3D magnetic field            */
    B_field_type_STS, /**< Spline-interpolated stellarator magnetic field   */
    B_field_type_TC   /**< Trivial Cartesian magnetic field                 */
} B_field_type;

/**
 * @brief Magnetic field offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized in
 * B_field_init_offload().
 *
 * The intended usage is that only single offload data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    B_field_type type;        /**< Magnetic field type wrapped by this struct */
    B_GS_offload_data BGS;    /**< GS field or NULL if not active             */
    B_2DS_offload_data B2DS;  /**< 2DS field or NULL if not active            */
    B_3DS_offload_data B3DS;  /**< 3DS field or NULL if not active            */
    B_STS_offload_data BSTS;  /**< STS field or NULL if not active            */
    B_TC_offload_data BTC;    /**< TC field or NULL if not active             */
    int offload_array_length; /**< Allocated offload array length             */
} B_field_offload_data;

/**
 * @brief Magnetic field simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from the B_field_offload_data in B_field_init().
 *
 * The intended usage is that only single B_field_data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    B_field_type type; /**< Magnetic field type wrapped by this struct */
    B_GS_data BGS;     /**< GS field or NULL if not active             */
    B_2DS_data B2DS;   /**< 2DS field or NULL if not active            */
    B_3DS_data B3DS;   /**< 3DS field or NULL if not active            */
    B_STS_data BSTS;   /**< STS field or NULL if not active            */
    B_TC_data BTC;     /**< TC field or NULL if not active             */
} B_field_data;

int B_field_init_offload(B_field_offload_data* offload_data,
                         real** offload_array);
void B_field_free_offload(B_field_offload_data* offload_data,
                          real** offload_array);

int B_field_init(
    B_field_data* Bdata, B_field_offload_data* offload_data,
    real* offload_array);
#ifndef GPU
#pragma omp declare simd uniform(Bdata)
#else
DECLARE_TARGET
#endif
a5err B_field_eval_psi(
    real* psi, real r, real phi, real z, real t, B_field_data* Bdata);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(Bdata)
#else
DECLARE_TARGET
#endif
a5err B_field_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, real t, B_field_data* Bdata);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(Bdata)
#else
DECLARE_TARGET
#endif
a5err B_field_eval_rho(real rho[2], real psi, B_field_data* Bdata);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(Bdata)
#else
DECLARE_TARGET
#endif
a5err B_field_eval_rho_drho(
    real rho_drho[4], real r, real phi, real z, B_field_data* Bdata);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(Bdata)
#else
DECLARE_TARGET
#endif
a5err B_field_eval_B(real B[3], real r, real phi, real z, real t,
                     B_field_data* Bdata);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(Bdata)
#else
DECLARE_TARGET
#endif
a5err B_field_eval_B_dB(
    real B_dB[15], real r, real phi, real z, real t, B_field_data* Bdata);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(Bdata)
#else
DECLARE_TARGET
#endif
a5err B_field_get_axis_rz(real rz[2], B_field_data* Bdata, real phi);
DECLARE_TARGET_END

#endif
