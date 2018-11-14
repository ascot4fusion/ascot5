/**
 * @file plasma.h
 * @brief Header file for plasma.c
 *
 * Contains a list declaring all plasma data, and declaration of
 * plasma_offload_data and plasma_data structs.
 */
#ifndef PLASMA_H
#define PLASMA_H

#include "ascot5.h"
#include "error.h"
#include "plasma/plasma_1D.h"
#include "plasma/plasma_1DS.h"

/**
 * @brief Plasma data types
 *
 * Plasma data types are used in the plassma interface (plasma.c) to direct
 * function calls to correct plasma data instances. Each plasma data  instance
 * must have a corresponding type.
 */
typedef enum plasma_type {
    plasma_type_1D, /**< Linear-interpolated 1D plasma data */
    plasma_type_1DS /**< Spline-interpolated 1D plasma data */
} plasma_type;


/**
 * @brief Plasma offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized in
 * plasma_init_offload().
 *
 * The intended usage is that only single offload data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    plasma_type
    type;                     /**< Plasma data type wrapped by this struct */
    plasma_1D_offload_data
    plasma_1D;                /**< 1D data or NULL if not active           */
    plasma_1DS_offload_data
    plasma_1DS;               /**< 1DS data or NULL if not active          */
    int offload_array_length; /**< Allocated offload array length          */
} plasma_offload_data;

/**
 * @brief Plasma simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from the plasma_offload_data in plasma_init().
 *
 * The intended usage is that only single plasma_data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    plasma_type type;           /**< Plasma data type wrapped by this struct */
    plasma_1D_data plasma_1D;   /**< 1D data or NULL if not active           */
    plasma_1DS_data plasma_1DS; /**< 1DS data or NULL if not active          */
} plasma_data;

int plasma_init_offload(plasma_offload_data* offload_data,
                        real** offload_array);
void plasma_free_offload(plasma_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
int plasma_init(plasma_data* plasma_data, plasma_offload_data* offload_data,
                real* offload_array);
#pragma omp declare simd uniform(pls_data)
real plasma_eval_temp(real rho, int species, plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data)
real plasma_eval_dens(real rho, int species, plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data)
a5err plasma_eval_densandtemp(
    real rho, plasma_data* pls_data, real* dens, real* temp);
#pragma omp declare simd uniform(pls_data)
int plasma_get_n_species(plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data)
real* plasma_get_species_mass(plasma_data* pls_data);
#pragma omp declare simd uniform(pls_data)
real* plasma_get_species_charge(plasma_data* pls_data);
#pragma omp end declare target

#endif
