/**
 * @author Joona Kontula joona.kontula@aalto.fi
 * @file E_1DS.h
 * @brief Header file for E_1DS.c
 *
 * Contains declaration of E_1DS_field_offload_data and E_1DS_field_data
 * structs.
 */
#ifndef E_1DS_H
#define E_1DS_H
#include "../offload_acc_omp.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"
#include "../B_field.h"

/**
 * @brief 1D spline electric field parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;                /**< Number of rho grid points           */
    real rho_min;             /**< Minimum rho value in the grid       */
    real rho_max;             /**< Maximum rho value in the grid       */
    int offload_array_length; /**< Number of elements in offload_array */
} E_1DS_offload_data;

/**
 * @brief 1D spline electric field parameters on the target
 */
typedef struct {
    interp1D_data dV;  /**< dV_drho 1D linear interpolation struct */
} E_1DS_data;

int E_1DS_init_offload(E_1DS_offload_data* offload_data, real** offload_array);
void E_1DS_free_offload(E_1DS_offload_data* offload_data, real** offload_array);

void E_1DS_init(E_1DS_data* Edata, E_1DS_offload_data* offload_data,
                real* offload_array);
#ifndef GPU
DECLARE_TARGET_SIMD_UNIFORM(Edata,Bdata)
#else
DECLARE_TARGET
#endif
a5err E_1DS_eval_E(real E[3], real r, real phi, real z, E_1DS_data* Edata,
                   B_field_data* Bdata);
DECLARE_TARGET_END
#endif
