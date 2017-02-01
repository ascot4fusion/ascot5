/**
 * @file B_field.h
 * @brief Header file for B_field.c
*/
#ifndef B_FIELD_H
#define B_FIELD_H

#include "B_GS.h"
#include "B_2D.h"
#include "B_3D.h"

typedef struct {
    int type;
    B_GS_offload_data BGS;
    B_2D_offload_data B2D;
    B_3D_offload_data B3D;
    int offload_array_length;
} B_field_offload_data;

typedef struct {
    int type;
    B_GS_data BGS;
    B_2D_data B2D;
    B_3D_data B3D;
} B_field_data;

void B_field_init_offload(B_field_offload_data* offload_data,
                          real** offload_array);
void B_field_free_offload(B_field_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
void B_field_init(B_field_data* Bdata, B_field_offload_data* offload_data,
                  real* offload_array);
#pragma omp declare simd uniform(Bdata)
void B_field_eval_B(real B[], real r, real phi, real z, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_field_eval_psi(real psi[], real r, real phi, real z,
                      B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_field_eval_rho(real rho[], real psi, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_field_eval_B_dB(real B_dB[], real r, real phi, real z,
                       B_field_data* Bdata);
real B_field_get_axis_r(B_field_data* Bdata);
real B_field_get_axis_z(B_field_data* Bdata);
#pragma omp end declare target   

#endif
