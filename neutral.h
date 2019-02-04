/**
 * @file neutral.h
 * @brief Header file for neutral.c
 *
 * Contains a list declaring all neutral_types, and declaration of
 * neutral_offload_data and neutral_data structs.
 */
#ifndef NEUTRAL_H
#define NEUTRAL_H

#include "ascot5.h"
#include "error.h"
#include "neutral/N0_3D.h"
#include "neutral/N0_ST.h"

/**
 * @brief Neutral data types
 *
 * Neutral data types are used in the neutral interface (neutral.c) to direct
 * function calls to correct neutral data instances. Each neutral data instance
 * must have a corresponding type.
 */
typedef enum neutral_type {
    neutral_type_3D, /**< Linearly-interpolated 3D neutral data          */
    neutral_type_ST  /**< Linearly-interpolated stellarator neutral data */
} neutral_type;

/**
 * @brief Neutral offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized in
 * neutral_init_offload().
 *
 * The intended usage is that only single offload data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    neutral_type type;        /**< Neutral data type wrapped by this struct */
    N0_3D_offload_data N03D;  /**< 3D field or NULL if not active           */
    N0_ST_offload_data N0ST;  /**< ST field or NULL if not active           */
    int offload_array_length; /**< Allocated offload array length           */
} neutral_offload_data;

/**
 * @brief Neutral simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from the neutral_offload_data in neutral_init().
 *
 * The intended usage is that only single neutral_data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    neutral_type type; /**< Neutral data type wrapped by this struct */
    N0_3D_data N03D;   /**< 3D field or NULL if not active           */
    N0_ST_data N0ST;   /**< ST field or NULL if not active           */
} neutral_data;

int neutral_init_offload(neutral_offload_data* offload_data,
                         real** offload_array);
void neutral_free_offload(neutral_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
int neutral_init(neutral_data* ndata, neutral_offload_data* offload_data,
                 real* offload_array);
#pragma omp declare simd uniform(ndata)
a5err neutral_eval_n0(real n0[], real r, real phi, real z, real t,
                      neutral_data* ndata);
#pragma omp end declare target
#endif
