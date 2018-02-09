/**
 * @file B_field.h
 * @brief Header file for B_field.c
*/
#ifndef B_FIELD_H
#define B_FIELD_H

#include "ascot5.h"
#include "error.h"
#include "Bfield/B_GS.h"
#include "Bfield/B_2D.h"
#include "Bfield/B_2DS.h"
#include "Bfield/B_3D.h"
#include "Bfield/B_3DS.h"
#include "Bfield/B_ST.h"
#include "Bfield/B_STS.h"
#include "Bfield/B_TC.h"

typedef enum B_field_type {
    B_field_type_GS, B_field_type_2D, B_field_type_2DS, B_field_type_3D, B_field_type_3DS, B_field_type_ST, B_field_type_STS, B_field_type_TC
} B_field_type;

typedef struct {
    B_field_type type;
    B_GS_offload_data BGS;
    B_2D_offload_data B2D;
    B_2DS_offload_data B2DS;
    B_3D_offload_data B3D;
    B_3DS_offload_data B3DS;
    B_ST_offload_data BST;
    B_STS_offload_data BSTS;
    B_TC_offload_data BTC;
    int offload_array_length;
} B_field_offload_data;

typedef struct {
    B_field_type type;
    B_GS_data BGS;
    B_2D_data B2D;
    B_2DS_data B2DS;
    B_3D_data B3D;
    B_3DS_data B3DS;
    B_ST_data BST;
    B_STS_data BSTS;
    B_TC_data BTC;
} B_field_data;

void B_field_init_offload(B_field_offload_data* offload_data,
                          real** offload_array);
void B_field_free_offload(B_field_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
int B_field_init(B_field_data* Bdata, B_field_offload_data* offload_data,
		 real* offload_array);
#pragma omp declare simd uniform(Bdata)
a5err B_field_eval_psi(real psi[], real r, real phi, real z,
		       B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_field_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z,
			    B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_field_eval_rho(real rho[], real psi, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_field_eval_rho_drho(real rho_drho[], real r, real phi, real z,
			    B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_field_eval_B(real B[], real r, real phi, real z, B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_field_eval_B_dB(real B_dB[], real r, real phi, real z,
			B_field_data* Bdata);
#pragma omp declare simd linear(i) uniform(psi, Bdata)
a5err B_field_eval_psi_SIMD(int i, real psi[NSIMD], real r, real phi, real z,
		       B_field_data* Bdata);
#pragma omp declare simd linear(i) uniform(rho, Bdata)
a5err B_field_eval_rho_SIMD(int i, real rho[NSIMD], real psi, B_field_data* Bdata);
#pragma omp declare simd linear(i) uniform(B_dB, Bdata)
a5err B_field_eval_B_dB_SIMD(int i, real B_dB[12][NSIMD], real r, real phi, real z,
			B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_field_get_axis_r(B_field_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_field_get_axis_z(B_field_data* Bdata);
#pragma omp end declare target   

#endif
