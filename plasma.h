/**
 * @file plasma.h
 * @brief Header file for plasma.c
*/
#ifndef PLASMA_H
#define PLASMA_H

#include "ascot5.h"
#include "error.h"
#include "plasma/plasma_1D.h"
#include "plasma/plasma_1DS.h"

/**
 * @brief All plasma input types.
 */
typedef enum plasma_type {
    plasma_type_1D, plasma_type_1DS
} plasma_type;

/**
 * @brief Plasma parameters that will be offloaded to target.
 */
typedef struct {
    plasma_type type;
    plasma_1D_offload_data plasma_1D;
    plasma_1DS_offload_data plasma_1DS;
    int offload_array_length;
} plasma_offload_data;

/**
 * @brief Plasma parameters on the target.
 */
typedef struct {
    plasma_type type;
    plasma_1D_data plasma_1D;
    plasma_1DS_data plasma_1DS;
} plasma_data;

void plasma_init_offload(plasma_offload_data* offload_data,
                          real** offload_array);
void plasma_free_offload(plasma_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
int plasma_init(plasma_data* plasma_data, plasma_offload_data* offload_data,
		 real* offload_array);
#pragma omp declare simd uniform(pls_data) simdlen(8)
real plasma_eval_temp(real rho, int species, plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data) simdlen(8)
real plasma_eval_dens(real rho, int species, plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data) simdlen(8)
a5err plasma_eval_densandtemp(real rho, plasma_data* pls_data, real* dens, real* temp);
#pragma omp declare simd uniform(pls_data) simdlen(8)
int plasma_get_n_species(plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data) simdlen(8)
real* plasma_get_species_mass(plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data) simdlen(8)
real* plasma_get_species_charge(plasma_data* pls_data);
#pragma omp end declare target   

#endif
