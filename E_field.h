/**
 * @file E_field.h
 * @brief Header file for E_field.c
*/
#ifndef E_FIELD_H
#define E_FIELD_H

#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "Efield/E_TC.h"
#include "Efield/E_1D.h"
#include "Efield/E_1DS.h"
#include "Efield/E_3D.h"

typedef enum E_field_type {
    E_field_type_TC, E_field_type_1D, E_field_type_1DS, E_field_type_3D
} E_field_type;

typedef struct {
    E_field_type type;
    E_TC_offload_data ETC;
    E_1D_offload_data E1D;
    E_1DS_offload_data E1DS;
    E_3D_offload_data E3D;
    int offload_array_length;
} E_field_offload_data;

typedef struct {
    E_field_type type;
    E_TC_data ETC;
    E_1D_data E1D;
    E_1DS_data E1DS;
    E_3D_data E3D;
} E_field_data;

void E_field_init_offload(E_field_offload_data* offload_data,
                          real** offload_array);
void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
int E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
		 real* offload_array);
#pragma omp declare simd uniform(Edata, Bdata)
a5err E_field_eval_E(real* E, real r, real phi, real z, E_field_data* Edata, B_field_data* Bdata);
#pragma omp declare simd uniform(Edata, Bdata)
a5err E_field_eval_E_SIMD(int i, real E[3][NSIMD], real r, real phi, real z, E_field_data* Edata, B_field_data* Bdata);
#pragma omp end declare target

#endif
