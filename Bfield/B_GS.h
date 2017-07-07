/**
 * @file B_GS.h
 * @brief Header file for B_GS.c
 */
#ifndef B_GS_H
#define B_GS_H
#include "ascot5.h"

/**
 * @brief Analytic magnetic field parameters that will be offloaded to target
 */
typedef struct {
    real psi0;                  /**< sqrt(psi) value at magnetic axis */
    real psi1;                  /**< sqrt(psi) value at separatrix */
    int offload_array_length;   /**< number of elements in offload_array */
} B_GS_offload_data;

/**
 * @brief Analytic magnetic field parameters on the target
 */
typedef struct {
    real psi0;          /**< sqrt(psi) value at magnetic axis */
    real psi1;          /**< sqrt(psi) value at separatrix */
    real *psi_coeff;    /**< Coefficients for the psi function components */
} B_GS_data;

void B_GS_init_offload(B_GS_offload_data* offload_data, real** offload_array);
void B_GS_free_offload(B_GS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_GS_init(B_GS_data* Bdata, B_GS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
void B_GS_eval_B(real B[], real r, real phi, real z, B_GS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_GS_eval_psi(real psi[], real r, real phi, real z, B_GS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_GS_eval_psi_dpsi(real psi[], real r, real phi, real z, B_GS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_GS_eval_rho(real rho[], real psi, B_GS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_GS_eval_rho_drho(real rho_drho[], real r, real phi, real z,
                        B_GS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_GS_eval_B_dB(real B_dB[], real r, real phi, real z, B_GS_data* Bdata);
real B_GS_get_axis_r(B_GS_data* Bdata);
real B_GS_get_axis_z(B_GS_data* Bdata);
#pragma omp end declare target

#endif
