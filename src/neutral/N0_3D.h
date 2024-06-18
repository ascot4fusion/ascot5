/**
 * @file N0_3D.h
 * @brief Header file for N0_3D.c
 */
#ifndef N0_3D_H
#define N0_3D_H
#include "../ascot5.h"
#include "../linint/linint.h" /* for 3D interpolation routines */

/**
 * @brief 3D neutral parameters on the host
 */
typedef struct {
    int n_r;        /**< Number of r grid points in the data                  */
    int n_z;        /**< Number of z grid points in the data                  */
    int n_phi;      /**< Number of phi grid points in the data                */
    real r_min;     /**< Minimum r coordinate in the grid in the data [m]     */
    real r_max;     /**< Maximum r coordinate in the grid in the data [m]     */
    real z_min;     /**< Minimum z coordinate in the grid in the data [m]     */
    real z_max;     /**< Maximum z coordinate in the grid in the data [m]     */
    real phi_min;   /**< Minimum phi coordinate in the grid in the data [rad] */
    real phi_max;   /**< Maximum phi coordinate in the grid in the data [rad] */
    int n_species;               /**< Number of neutral species               */
    int anum[MAX_SPECIES];       /**< Neutral species mass number []          */
    int znum[MAX_SPECIES];       /**< Neutral species charge number []        */
    int maxwellian[MAX_SPECIES]; /**< Whether species distribution is
                                    Maxwellian or monoenergetic               */
    int offload_array_length;    /**< Number of elements in offload_array     */
} N0_3D_offload_data;

/**
 * @brief 3D neutral parameters on the target
 */
typedef struct {
    int n_species;            /**< Number of neutral species                  */
    int anum[MAX_SPECIES];    /**< Neutral species mass number []             */
    int znum[MAX_SPECIES];    /**< Neutral species charge number []           */
    int maxwellian[MAX_SPECIES];    /**< Whether species distribution is
                                       Maxwellian or monoenergetic            */
    linint3D_data n0[MAX_SPECIES];  /**< Pointer to start of neutral density
                                       interpolation data struct array        */
    linint3D_data t0[MAX_SPECIES];  /**< Pointer to start of neutral temperature
                                       interpolation data struct array        */
} N0_3D_data;

int N0_3D_init_offload(N0_3D_offload_data* offload_data, real** offload_array);
void N0_3D_free_offload(N0_3D_offload_data* offload_data, real** offload_array);

void N0_3D_init(N0_3D_data* ndata, N0_3D_offload_data* offload_data,
                real* offload_array);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_3D_eval_n0(real* n0, real r, real phi, real z, N0_3D_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_3D_eval_t0(real* t0, real r, real phi, real z, N0_3D_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
int N0_3D_get_n_species(N0_3D_data* ndata);
#endif
