/**
 * @file N0_1D.h
 * @brief Header file for N0_1D.c
 */
#ifndef N0_1D_H
#define N0_1D_H
#include "../ascot5.h"
#include "../linint/linint.h" /* for 1D interpolation routines */

/**
 * @brief 1D neutral parameters on the host
 */
typedef struct {
    int n_rho;      /**< Number of r grid points in the data                  */
    real rho_min;   /**< Minimum r coordinate in the grid in the data [m]     */
    real rho_max;   /**< Maximum r coordinate in the grid in the data [m]     */
    int n_species;               /**< Number of neutral species               */
    int anum[MAX_SPECIES];       /**< Neutral species mass number []          */
    int znum[MAX_SPECIES];       /**< Neutral species charge number []        */
    int maxwellian[MAX_SPECIES]; /**< Whether species distribution is
                                    Maxwellian or monoenergetic               */
    int offload_array_length;    /**< Number of elements in offload_array     */
} N0_1D_offload_data;

/**
 * @brief 1D neutral parameters on the target
 */
typedef struct {
    int n_species;            /**< Number of neutral species                  */
    int anum[MAX_SPECIES];    /**< Neutral species mass number []             */
    int znum[MAX_SPECIES];    /**< Neutral species charge number []           */
    int maxwellian[MAX_SPECIES];    /**< Whether species distribution is
                                       Maxwellian or monoenergetic            */
    linint1D_data n0[MAX_SPECIES];  /**< Pointer to start of neutral density
                                       interpolation data struct array        */
    linint1D_data t0[MAX_SPECIES];  /**< Pointer to start of neutral temperature
                                       interpolation data struct array        */
} N0_1D_data;

int N0_1D_init_offload(N0_1D_offload_data* offload_data, real** offload_array);
void N0_1D_free_offload(N0_1D_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void N0_1D_init(N0_1D_data* ndata, N0_1D_offload_data* offload_data,
                real* offload_array);
#pragma omp declare simd uniform(ndata)
a5err N0_1D_eval_n0(real* n0, real rho, int species,
                    N0_1D_data* ndata);
#pragma omp declare simd uniform(ndata)
a5err N0_1D_eval_t0(real* t0, real rho, int species,
                    N0_1D_data* ndata);
#pragma omp end declare target
#endif
