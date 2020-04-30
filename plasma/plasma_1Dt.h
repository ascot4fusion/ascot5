/**
 * @file plasma_1Dt.h
 * @brief Header file for plasma_1Dt.c
 */
#ifndef PLASMA_1DT_H
#define PLASMA_1DT_H
#include "../ascot5.h"
#include "../error.h"

/**
 * @brief 1D plasma parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data    */
    int n_time;                 /**< number of time points               */
    int n_species;              /**< number of plasma species including
                                     electrons                           */
    real mass[MAX_SPECIES];     /**< plasma species masses [kg]          */
    real charge[MAX_SPECIES];   /**< plasma species charges [C]          */
    int anum[MAX_SPECIES];      /**< ion species atomic number           */
    int znum[MAX_SPECIES];      /**< ion species charge number           */
    int offload_array_length;   /**< number of elements in offload_array */
} plasma_1Dt_offload_data;

/**
 * @brief 1D plasma parameters on the target
 */
typedef struct {
    int n_rho;                /**< number of rho values in the data    */
    int n_time;               /**< number of time points               */
    int n_species;            /**< number of plasma species including
                                   electrons                           */
    real mass[MAX_SPECIES];   /**< plasma species masses [kg]          */
    real charge[MAX_SPECIES]; /**< plasma species charges [C]          */
    int anum[MAX_SPECIES];    /**< ion species atomic number           */
    int znum[MAX_SPECIES];    /**< ion species charge number           */
    real* rho;                /**< pointer to start of rho values in
                                   offload_array                       */
    real* time;               /**< pointer to start of time values     */
    real* temp;               /**< pointer to start of temperatures    */
    real* dens;               /**< pointer to start of densities       */
} plasma_1Dt_data;

int plasma_1Dt_init_offload(plasma_1Dt_offload_data* offload_data,
                           real** offload_array);
void plasma_1Dt_free_offload(plasma_1Dt_offload_data* offload_data,
                            real** offload_array);

#pragma omp declare target
void plasma_1Dt_init(plasma_1Dt_data* pls_data,
                    plasma_1Dt_offload_data* offload_data,
                    real* offload_array);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1Dt_eval_temp(real* dens, real rho, int species,
                          plasma_1Dt_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1Dt_eval_dens(real* temp, real rho, int species,
                          plasma_1Dt_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1Dt_eval_densandtemp(real* dens, real* temp, real rho,
                                 plasma_1Dt_data* pls_data);
#pragma omp declare simd uniform(pls_data)
int plasma_1Dt_get_n_species(plasma_1Dt_data* pls_data);
#pragma omp declare simd uniform(pls_data)
const real* plasma_1Dt_get_species_mass(plasma_1Dt_data* pls_data);
#pragma omp declare simd uniform(pls_data)
const real* plasma_1Dt_get_species_charge(plasma_1Dt_data* pls_data);
#pragma omp end declare target

#endif
