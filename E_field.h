/**
 * @file E_field.h
 * @brief Header file for E_field.c
*/
#ifndef E_FIELD_H
#define E_FIELD_H

#include "E_1D.h"

typedef enum E_field_type {
    E_field_type_1D
} E_field_type;

typedef struct {
    E_field_type type;
    E_1D_offload_data E1D;
    int offload_array_length;
} E_field_offload_data;

typedef struct {
    E_field_type type;
    E_1D_data E1D;
} E_field_data;

void E_field_init_offload(E_field_offload_data* offload_data,
                          real** offload_array);
void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
void E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
                  real* offload_array);
#pragma omp declare simd uniform(Edata) 
void E_field_eval_E(real* E, real* rho_drho, E_field_data* Edata);
#pragma omp end declare target

#endif
