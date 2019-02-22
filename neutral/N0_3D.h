/**
 * @file N0_3D.h
 * @brief Header file for N0_3D.c
 */
#ifndef N0_3D_H
#define N0_3D_H
#include "../ascot5.h"
#include "../linint/linint.h" /* for 3D interpolation routines */

/**
 * @brief 3D neutral density parameters on the host
 */
typedef struct {
    int n_r;        /**< number of r grid points in the data                  */
    int n_z;        /**< number of z grid points in the data                  */
    int n_phi;      /**< number of phi grid points in the data                */
    real r_min;     /**< minimum r coordinate in the grid in the data         */
    real r_max;     /**< maximum r coordinate in the grid in the data         */
    real r_grid;    /**< r grid interval (r_max-r_min)/(n_r-1) in the data    */
    real z_min;     /**< minimum z coordinate in the grid in the data         */
    real z_max;     /**< maximum z coordinate in the grid in the data         */
    real z_grid;    /**< z grid interval (z_max-z_min)/(n_z-1) in the data    */
    real phi_min;   /**< minimum phi coordinate in the grid in the data       */
    real phi_max;   /**< maximum phi coordinate in the grid in the data       */
    real phi_grid;  /**< phi grid interval 2pi/(n_phi-1) in the data          */
    int n_species;              /**< number of neutral species                */
    real mass[MAX_SPECIES];     /**< neutral species masses [kg]              */
    real anum[MAX_SPECIES];     /**< neutral species mass number []           */
    real znum[MAX_SPECIES];     /**< neutral species charge number []         */
    int offload_array_length;   /**< number of elements in offload_array      */
} N0_3D_offload_data;

/**
 * @brief 3D neutral density parameters on the target
 */
typedef struct {
    int n_species;                 /**< number of neutral species             */
    real mass[MAX_SPECIES];        /**< neutral species masses [kg]           */
    real anum[MAX_SPECIES];        /**< neutral species mass number []        */
    real znum[MAX_SPECIES];        /**< neutral species charge number []      */
    linint3D_data n0[MAX_SPECIES]; /**< pointer to start of neutral density
                                      interpolation data struct array         */
} N0_3D_data;

int N0_3D_init_offload(N0_3D_offload_data* offload_data, real** offload_array);
void N0_3D_free_offload(N0_3D_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void N0_3D_init(N0_3D_data* ndata, N0_3D_offload_data* offload_data,
                real* offload_array);
#pragma omp declare simd uniform(ndata)
a5err N0_3D_eval_n0(real* n0, real r, real phi, real z, int species,
                    N0_3D_data* ndata);
#pragma omp end declare target
#endif
