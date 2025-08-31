/**
 * @file B_field.h
 * @brief Header file for B_field.c
 *
 * Contains a list declaring all B_field_types, and declaration of
 * B_field_offload_data and B_field_data structs.
 */
#ifndef B_FIELD_H
#define B_FIELD_H

#include "ascot5.h"
#include "offload.h"
#include "error.h"
#include "Bfield/B_GS.h"
#include "Bfield/B_2DS.h"
#include "Bfield/B_3DS.h"
#include "Bfield/B_STS.h"
#include "Bfield/B_TC.h"

/**
 * @brief Magnetic field types
 *
 * Magnetic field types are used in the magnetic field interface to direct
 * function calls to correct magnetic field instances. Each magnetic field
 * instance must have a corresponding type.
 */
typedef enum B_field_type {
    B_field_type_TC,  /**< Trivial Cartesian magnetic field                 */
    B_field_type_GS,  /**< Analytic magnetic field                          */
    B_field_type_2DS, /**< Spline-interpolated axisymmetric  magnetic field */
    B_field_type_3DS, /**< Spline-interpolated 3D magnetic field            */
    B_field_type_STS, /**< Spline-interpolated stellarator magnetic field   */
} B_field_type;

/**
 * @brief Magnetic field simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct {
    B_TC_data* BTC;    /**< TC field or NULL if not active             */
    B_GS_data* BGS;    /**< GS field or NULL if not active             */
    B_2DS_data* B2DS;  /**< 2DS field or NULL if not active            */
    B_3DS_data* B3DS;  /**< 3DS field or NULL if not active            */
    B_STS_data* BSTS;  /**< STS field or NULL if not active            */
    B_field_type type; /**< Magnetic field type wrapped by this struct */
} B_field_data;

void B_field_free(B_field_data* data);
void B_field_offload(B_field_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_field_eval_psi(
    real* psi, real r, real phi, real z, real t, B_field_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_field_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, real t, B_field_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_field_eval_rho(real rho[2], real psi, B_field_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_field_eval_rho_drho(
    real rho_drho[4], real r, real phi, real z, B_field_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_field_eval_B(real B[3], real r, real phi, real z, real t,
                     B_field_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_field_eval_B_dB(
    real B_dB[15], real r, real phi, real z, real t, B_field_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(Bdata)
a5err B_field_get_axis_rz(real rz[2], B_field_data* Bdata, real phi);
DECLARE_TARGET_END

#endif
