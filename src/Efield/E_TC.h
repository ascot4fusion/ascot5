/**
 * @file E_TC.h
 * @brief Header file for E_TC.c
 *
 * Contains declaration of E_TC_field_offload_data and E_TC_field_data structs.
 */
#ifndef E_TC_H
#define E_TC_H

#include "../offload_acc_omp.h"
#include "../ascot5.h"
#include "../error.h"
#include "../B_field.h"

/**
 * @brief Trivial Cartesian electric field offload data
 */
typedef struct {
    real Exyz[3]; /**< Pointer to array holding constant [E_x, E_y, E_z]
                     values [V/m]                                           */
    int offload_array_length; /**< Offload array length, zero        always */
} E_TC_offload_data;

/**
 * @brief Trivial Cartesian electric field simulation data
 */
typedef struct {
    real* Exyz; /**< Pointer to array holding constant [E_x, E_y, E_z]
                     values [V/m]                                             */
} E_TC_data;

int E_TC_init_offload(E_TC_offload_data* offload_data,
                      real** offload_array);
void E_TC_free_offload(E_TC_offload_data* offload_data,
                       real** offload_array);


void E_TC_init(E_TC_data* Edata, E_TC_offload_data* offload_data,
                  real* offload_array);
#ifndef GPU
#pragma omp declare simd uniform(Edata,Bdata)
#else
DECLARE_TARGET
#endif
a5err E_TC_eval_E(real E[3], real r, real phi, real z, E_TC_data* Edata,
                  B_field_data* Bdata);
DECLARE_TARGET_END
#endif
