/**
 * @file plasma_1D.h
 * @brief Header file for plasma_1D.c
 */
#ifndef PLASMA_1D_H
#define PLASMA_1D_H
#include "../ascot5.h"
#include "../error.h"

/**
 * @brief 1D plasma parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data    */
    int n_species;              /**< number of plasma species including
                                     electrons                           */
    real mass[MAX_SPECIES];     /**< plasma species masses [kg]          */
    real charge[MAX_SPECIES];   /**< plasma species charges [C]          */
    int anum[MAX_SPECIES];      /**< ion species atomic number           */
    int znum[MAX_SPECIES];      /**< ion species charge number           */
    int offload_array_length;   /**< number of elements in offload_array */
} plasma_1D_offload_data;

/**
 * @brief 1D plasma parameters on the target
 */
typedef struct {
    int n_rho;                /**< number of rho values in the data      */
    int n_species;            /**< number of plasma species including
                                   electrons                             */
    real mass[MAX_SPECIES];   /**< plasma species masses [kg]            */
    real charge[MAX_SPECIES]; /**< plasma species charges [C]            */
    int anum[MAX_SPECIES];    /**< ion species atomic number             */
    int znum[MAX_SPECIES];    /**< ion species charge number             */
    real* rho;                /**< pointer to start of rho values in
                                   offload_array                         */
    real* temp;               /**< pointer to start of temperatures      */
    real* dens;               /**< pointer to start of densities         */
    real* vtor;               /**< pointer to start of toroidal rotation */
} plasma_1D_data;

int plasma_1D_init_offload(plasma_1D_offload_data* offload_data,
                           real** offload_array);
void plasma_1D_free_offload(plasma_1D_offload_data* offload_data,
                            real** offload_array);

#pragma omp declare target
void plasma_1D_init(plasma_1D_data* pls_data,
                    plasma_1D_offload_data* offload_data,
                    real* offload_array);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1D_eval_temp(real* dens, real rho, int species,
                          plasma_1D_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1D_eval_dens(real* temp, real rho, int species,
                          plasma_1D_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1D_eval_densandtemp(real* dens, real* temp, real rho,
                                 plasma_1D_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1D_eval_flow(
    real* vflow, real rho, real r, plasma_1D_data* pls_data);
#pragma omp end declare target

#endif
