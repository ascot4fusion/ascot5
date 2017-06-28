/**
 * @file diag.h
 * @brief Header file for diag.c
 */
#ifndef DIAG_H
#define DIAG_H
#include "ascot5.h"
#include "particle.h"
#include "distributions.h"
#include "diag_orb.h"

typedef struct{
    int orb_collect;
    int debug_collect;
    int dist4D_collect;
    diag_orb_offload_data orbits;
    dist_rzvv_offload_data dist4D;

    int offload_dist4D_index;
    int offload_array_length; /**< number of elements in offload_array */

} diag_offload_data;

typedef struct{
    int diag_orb_collect;
    int diag_debug_collect;
    int diag_dist4D_collect;
    diag_orb_data orbits;
    dist_rzvv_data dist4D;

    int offload_dist4D_index;

} diag_data;

void diag_init_offload(diag_offload_data* data, real** offload_array);

void diag_free_offload(diag_offload_data* data, real** offload_array);

void diag_sum(diag_data* d, real* array1, real* array2);

#pragma omp declare target
void diag_init(diag_data* data, diag_offload_data* offload_data, real* offload_array);

void diag_update_gc(diag_data* d, particle_simd_gc* p_f, particle_simd_gc* p_i);

void diag_update_fo(diag_data* d, particle_simd_fo* p_f, particle_simd_fo* p_i);

void diag_update_ml(diag_data* d, particle_simd_ml* p_f, particle_simd_ml* p_i);

void diag_clean(diag_data* d);
#pragma omp end declare target

#endif

