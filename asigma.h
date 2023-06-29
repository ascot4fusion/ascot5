/**
 * @file asigma.h
 * @brief Header file for asigma.c
 *
 * Contains a list declaring all atomic reaction data types, and declaration
 * of asigma_offload_data and asigma_data structs.
 */
#ifndef ASIGMA_H
#define ASIGMA_H

#include "ascot5.h"
#include "error.h"

/**
 * @brief Atomic reaction data types
 *
 * Atomic reaction data types are used in the atomic reaction data interface
 * (asigma.c) to direct function calls to correct atomic reaction data
 * instances. Each atomic reaction data instance must have a
 * corresponding type.
 */
typedef enum asigma_type {
    asigma_type_loc, /**< Atomic reaction data from local files */
} asigma_type;

/**
 * @brief Atomic reaction offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized
 * in asigma_init_offload().
 */
typedef struct {
    asigma_type type;   /**< Atomic reaction data type wrapped by this struct */
    int offload_array_length; /**< Allocated offload array length */
} asigma_offload_data;

/**
 * @brief Atomic reaction simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from input data in asigma_init().
 */
typedef struct {
    asigma_type type; /**< Atomic reaction data type wrapped by this struct */
} asigma_data;

int asigma_init_offload(asigma_offload_data* offload_data,
                        real** offload_array);
void asigma_free_offload(asigma_offload_data* offload_data,
                         real** offload_array);

#pragma omp declare target
int asigma_init(asigma_data* asigma_data, asigma_offload_data* offload_data,
                real* offload_array);
#pragma omp declare simd uniform(asigmadata)
a5err asigma_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, int reac_type,
    asigma_data* asigmadata, real E_coll_per_amu);
#pragma omp declare simd uniform(asigmadata)
a5err asigma_eval_sigmav(
    real* sigmav, int z_1, int a_1, int z_2, int a_2, int reac_type,
    asigma_data* asigmadata, real E, real T_e, real* T_i, real T_0,
    real n_e, real* n_i, int i_spec);
#pragma omp end declare target

#endif
