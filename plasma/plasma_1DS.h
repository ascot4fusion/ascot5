/**
 * @file plasma_1DS.h
 * @brief Header file for plasma_1DS.c
 */
#ifndef PLASMA_1DS_H
#define PLASMA_1DS_H
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp1Dcomp.h"

/**
 * @brief 1D plasma parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data */
    int n_species;              /**< number of plasma species; first is
                                     electrons, then ions */
    real mass[MAX_SPECIES];     /**< plasma species masses (kg) */
    real charge[MAX_SPECIES];   /**< plasma species charges (C) */
    int offload_array_length;   /**< number of elements in offload_array */
} plasma_1DS_offload_data;

/**
 * @brief 1D plasma parameters on the target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data */
    real rho_min;                /**< number of rho values in the data */
    real rho_max;                /**< number of rho values in the data */
    real rho_grid;                /**< number of rho values in the data */
    int n_species;              /**< number of plasma species; first is
                                     electrons, then ions */
    real mass[MAX_SPECIES];     /**< plasma species masses (kg) */
    real charge[MAX_SPECIES];   /**< plasma species charges (C) */
    interp1D_data temp[2];         /**< pointer to start of temperature interpolation data structs */
    interp1D_data* dens;                 /**< pointer to start of densities */
} plasma_1DS_data;

void plasma_1DS_free(plasma_1DS_data* pls_data);

void plasma_1DS_init_offload(plasma_1DS_offload_data* offload_data,
                             real** offload_array);

void plasma_1DS_free_offload(plasma_1DS_offload_data* offload_data,
                             real** offload_array);

#pragma omp declare target
a5err plasma_1DS_init(plasma_1DS_data* pls_data,
                      plasma_1DS_offload_data* offload_data,
                      real* offload_array);

#pragma omp declare simd uniform(pls_data)
a5err plasma_1DS_eval_temp(real temp[], real rho, int species, plasma_1DS_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1DS_eval_dens(real dens[], real rho, int species, plasma_1DS_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_1DS_eval_densandtemp(real* dens, real* temp, real rho, plasma_1DS_data* pls_data);
#pragma omp end declare target

#endif
