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
#include "plasma/plasma_1Dt.h"
#include "plasma/plasma_1DS.h"

/**
 * @brief Plasma data types
 */
typedef enum plasma_type {
    plasma_type_1D,  /**< Linear-interpolated 1D plasma data                */
    plasma_type_1Dt, /**< Linear-interpolated time-dependent 1D plasma data */
    plasma_type_1DS  /**< Spline-interpolated 1D plasma data                */
} plasma_type;

/**
 * @brief Plasma simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct {
    plasma_type type;           /**< Plasma data type wrapped by this struct */
    plasma_1D_data plasma_1D;   /**< 1D data or NULL if not active           */
    plasma_1Dt_data plasma_1Dt; /**< 1D data or NULL if not active           */
    plasma_1DS_data plasma_1DS; /**< 1DS data or NULL if not active          */
} plasma_data;

void plasma_free(plasma_data* data);
void plasma_offload(plasma_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_eval_temp(real* temp, real rho, real r, real phi, real z, real t,
                       int species, plasma_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_eval_dens(real* dens, real rho, real r, real phi, real z, real t,
                       int species, plasma_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_eval_densandtemp(real* dens, real* temp, real rho,
                              real r, real phi, real z, real t,
                              plasma_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
a5err plasma_eval_flow(real* vflow, real rho, real r, real phi, real z, real t,
                       plasma_data* pls_data);
DECLARE_TARGET_END
int plasma_get_n_species(plasma_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
const real* plasma_get_species_mass(plasma_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
const real* plasma_get_species_charge(plasma_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
const int* plasma_get_species_znum(plasma_data* pls_data);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(pls_data)
const int* plasma_get_species_anum(plasma_data* pls_data);
DECLARE_TARGET_END

#endif
