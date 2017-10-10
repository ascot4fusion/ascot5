/**
 * @file E_TC.h
 * @brief Header file for E_TC.c
*/
#ifndef E_TC_H
#define E_TC_H

#include "../ascot5.h"
#include "../B_field.h"

typedef struct {
    int offload_array_length;
} E_TC_offload_data;

typedef struct {
    real* Exyz;
} E_TC_data;

void E_TC_init_offload(E_TC_offload_data* offload_data,
			     real** offload_array);
void E_TC_free_offload(E_TC_offload_data* offload_data,
			     real** offload_array);

#pragma omp declare target
void E_TC_init(E_TC_data* Edata, E_TC_offload_data* offload_data,
                  real* offload_array);
#pragma omp declare simd uniform(Edata,Bdata) simdlen(8)
void E_TC_eval_E(real* E, real r, real phi, real z, E_TC_data* Edata, B_field_data* Bdata);
#pragma omp declare simd uniform(Edata,Bdata) simdlen(8)
void E_TC_eval_E_SIMD(int i, real* E, real r, real phi, real z, E_TC_data* Edata, B_field_data* Bdata);
#pragma omp end declare target

#endif
