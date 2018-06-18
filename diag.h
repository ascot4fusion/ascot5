/**
 * @file diag.h
 * @brief Header file for diag.c
 */
#ifndef DIAG_H
#define DIAG_H
#include "ascot5.h"
#include "particle.h"
#include "dist_5D.h"
#include "dist_6D.h"
#include "dist_rho5D.h"
#include "dist_rho6D.h"
#include "diag_orb.h"

typedef struct {
    int orb_collect;
    int debug_collect;
    int dist5D_collect;
    int dist6D_collect;
    int distrho5D_collect;
	int distrho6D_collect;

    diag_orb_offload_data orbits;
    dist_5D_offload_data dist5D;
    dist_6D_offload_data dist6D;
    dist_rho5D_offload_data distrho5D;
    dist_rho6D_offload_data distrho6D;

    int offload_dist5D_index;
    int offload_dist6D_index;
	int offload_distrho5D_index;
	int offload_distrho6D_index;
    int offload_array_length; /**< number of elements in offload_array */
} diag_offload_data;

typedef struct {
    int diag_orb_collect;
    int diag_debug_collect;
    int diag_dist5D_collect;
    int diag_dist6D_collect;
    int diag_distrho5D_collect;
    int diag_distrho6D_collect;

    diag_orb_data orbits;
    dist_5D_data dist5D;
    dist_6D_data dist6D;
    dist_rho5D_data distrho5D;
    dist_rho6D_data distrho6D;
    
    int offload_dist5D_index;
    int offload_dist6D_index;
    int offload_distrho5D_index;
    int offload_distrho6D_index;
} diag_data;

/** @brief Struct for storing particle specific data needed exclusively for diagnostics
 *  In principle, this could be stored in diag_data struct bu we need multiple instances,
 *  one for each thread. */
typedef struct {
    integer particleId[NSIMD];
    real prevWriteTime[NSIMD];
    int nextN[NSIMD];
    diag_orb_dat** Nlist;
} diag_storage;

void diag_init_offload(diag_offload_data* data, real** offload_array);

void diag_free_offload(diag_offload_data* data, real** offload_array);

void diag_sum(diag_data* d, real* array1, real* array2);

#pragma omp declare target
void diag_init(diag_data* data, diag_offload_data* offload_data, real* offload_array);

void diag_update_fo(diag_data* d, diag_storage* ds, particle_simd_fo* p_f, particle_simd_fo* p_i);

void diag_update_gc(diag_data* d, diag_storage* ds, particle_simd_gc* p_f, particle_simd_gc* p_i);

void diag_update_ml(diag_data* d, diag_storage* ds, particle_simd_ml* p_f, particle_simd_ml* p_i);

void diag_storage_aquire(diag_data* data, diag_storage** ds);

void diag_storage_discard(diag_storage* ds);

void diag_clean(diag_data* d);
#pragma omp end declare target

#endif
