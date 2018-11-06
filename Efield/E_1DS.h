/**
 * @author Joona Kontula joona.kontula@aalto.fi
 * @file E_1DS.h
 * @brief Header file for E_1DS.c
 */
#ifndef E_1DS_H
#define E_1DS_H
#include "../ascot5.h"
#include "../spline/interp1Dcomp.h"
#include "../B_field.h"

/**
 * @brief 1D spline electric field parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data */
    real rho_min;               /**< minimum rho value in the data */
    real rho_max;               /**< maximum rho value in the data */
    int offload_array_length;   /**< number of elements in offload_array */
} E_1DS_offload_data;

/**
 * @brief 1D spline electric field parameters on the target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data */
    real rho_min;               /**< minimum rho value in the data */
    real rho_max;               /**< maximum rho value in the data */
    real rho_grid;               /**< maximum rho value in the data */
    interp1D_data dV;           /**< spline representation of potential derivative values */
} E_1DS_data;

void E_1DS_free_offload(E_1DS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
a5err E_1DS_init(E_1DS_data* Edata, E_1DS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Edata,Bdata)  
a5err E_1DS_eval_E(real E[], real r, real phi, real z, E_1DS_data* Edata, B_field_data* Bdata);
#pragma omp end declare target
#endif
