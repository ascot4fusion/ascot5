/**
 * @file plasma_1DS.h
 * @brief Header file for plasma_1DS.c
 */
#ifndef PLASMA_1DS_H
#define PLASMA_1DS_H
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 1D spline plasma parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;                  /**< number of rho values in the data    */
    real rho_min;               /**< minimum rho value in the grid       */
    real rho_max;               /**< maximum rho value in the grid       */
    int n_species;              /**< number of plasma species including
                                     electrons                           */
    real mass[MAX_SPECIES];     /**< plasma species masses (kg)          */
    real charge[MAX_SPECIES];   /**< plasma species charges (C)          */
    int anum[MAX_SPECIES];      /**< ion species atomic number           */
    int znum[MAX_SPECIES];      /**< ion species charge number           */
    int offload_array_length;   /**< number of elements in offload_array */
} plasma_1DS_offload_data;

/**
 * @brief 1D spline plasma parameters on the target
 */
typedef struct {
    int n_species;              /**< number of plasma species including
                                     electrons                                */
    real mass[MAX_SPECIES];     /**< plasma species masses (kg)               */
    real charge[MAX_SPECIES];   /**< plasma species charges (C)               */
    int anum[MAX_SPECIES];      /**< ion species atomic number                */
    int znum[MAX_SPECIES];      /**< ion species charge number                */
    interp1D_data temp[2];      /**< electron and ion temperature
                                     interpolation structs                    */
    interp1D_data dens[MAX_SPECIES]; /**< electron and every ion species
                                          density interpolation structs       */
} plasma_1DS_data;

int plasma_1DS_init_offload(plasma_1DS_offload_data* offload_data,
                            real** offload_array);

void plasma_1DS_free_offload(plasma_1DS_offload_data* offload_data,
                             real** offload_array);

void plasma_1DS_init(plasma_1DS_data* pls_data,
                     plasma_1DS_offload_data* offload_data,
                     real* offload_array);
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
a5err plasma_1DS_eval_temp(real* temp, real rho, int species,
                           plasma_1DS_data* pls_data);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
a5err plasma_1DS_eval_dens(real* dens, real rho, int species,
                           plasma_1DS_data* pls_data);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
a5err plasma_1DS_eval_densandtemp(real* dens, real* temp, real rho,
                                  plasma_1DS_data* pls_data);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
int plasma_1DS_get_n_species(plasma_1DS_data* pls_data);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
const real* plasma_1DS_get_species_mass(plasma_1DS_data* pls_data);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
const real* plasma_1DS_get_species_charge(plasma_1DS_data* pls_data);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
const int* plasma_1DS_get_species_znum(plasma_1DS_data* pls_data);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(pls_data)
#else
DECLARE_TARGET
#endif
const int* plasma_1DS_get_species_anum(plasma_1DS_data* pls_data);
DECLARE_TARGET_END

#endif
