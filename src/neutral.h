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
#include "neutral/N0_1D.h"
#include "neutral/N0_3D.h"

/**
 * @brief Neutral data types
 */
typedef enum neutral_type {
    neutral_type_1D, /**< Linearly-interpolated 1D neutral data          */
    neutral_type_3D, /**< Linearly-interpolated 3D neutral data          */
} neutral_type;

/**
 * @brief Neutral simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct {
    N0_1D_data* N01D;   /**< 1D field or NULL if not active           */
    N0_3D_data* N03D;   /**< 3D field or NULL if not active           */
    neutral_type type;  /**< Neutral data type wrapped by this struct */
} neutral_data;

void neutral_free(neutral_data* data);
void neutral_offload(neutral_data* data);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err neutral_eval_n0(real* n0, real rho, real r, real phi, real z, real t,
                      neutral_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err neutral_eval_t0(real* t0, real rho, real r, real phi, real z, real t,
                      neutral_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
int neutral_get_n_species(neutral_data* ndata);
#endif
