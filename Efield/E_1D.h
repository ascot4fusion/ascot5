/**
 * @author Joona Kontula joona.kontula@aalto.fi
 * @file E_1D.h
 * @brief Header file for E_1D.c
 */
#ifndef E_1D_H
#define E_1D_H
#include "../ascot5.h"
#include "../B_field.h"

/**
 * @brief 1D electric field parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data */
    real rho_min;               /**< minimum rho value in the data */
    real rho_max;               /**< maximum rho value in the data */
    int offload_array_length;   /**< number of elements in offload_array */
} E_1D_offload_data;

/**
 * @brief 1D electric field parameters on the target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data */
    real rho_min;               /**< minimum rho value in the data */
    real rho_max;               /**< maximum rho value in the data */
    real* rho;                  /**< pointer to start of rho values in 
                                   offload_array */
    real* dV;                   /**< pointer to start of potential derivative values */
} E_1D_data;

void E_1D_init_offload(E_1D_offload_data* offload_data, real** offload_array);
void E_1D_free_offload(E_1D_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void E_1D_init(E_1D_data* Edata, E_1D_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Edata)
void E_1D_eval_E(real E[], real r, real phi, real z, E_1D_data* Edata, B_field_data* Bdata);
#pragma omp end declare target   
#endif
