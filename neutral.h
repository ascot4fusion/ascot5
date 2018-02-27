/**
 * @file neutral.h
 * @brief Header file for neutral.c
*/
#ifndef NEUTRAL_H
#define NEUTRAL_H

#include "ascot5.h"
#include "error.h"
#include "neutral/N0_3D.h"
#include "neutral/N0_ST.h"

typedef enum neutral_type {
    neutral_type_3D, neutral_type_ST
} neutral_type;

typedef struct {
    neutral_type type;
    N0_3D_offload_data N03D;
    N0_ST_offload_data N0ST;
    int offload_array_length;
} neutral_offload_data;

typedef struct {
    neutral_type type;
    N0_3D_data N03D;
    N0_ST_data N0ST;
} neutral_data;

void neutral_init_offload(neutral_offload_data* offload_data,
                          real** offload_array);
void neutral_free_offload(neutral_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
int neutral_init(neutral_data* ndata, neutral_offload_data* offload_data,
		 real* offload_array);
#pragma omp declare simd uniform(ndata)
a5err neutral_eval_n0(real n0[], real r, real phi, real z,
		       neutral_data* ndata);
#pragma omp declare simd linear(i) uniform(n0, ndata)
a5err neutral_eval_n0_SIMD(int i, real n0[NSIMD], real r, real phi, real z,
		       neutral_data* ndata);
#pragma omp end declare target   
#endif
