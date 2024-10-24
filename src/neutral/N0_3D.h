/**
 * @file N0_3D.h
 * @brief Header file for N0_3D.c
 */
#ifndef N0_3D_H
#define N0_3D_H
#include "../ascot5.h"
#include "../offload.h"
#include "../linint/linint.h" /* for 3D interpolation routines */

/**
 * @brief 3D neutral parameters on the target
 */
typedef struct {
    int n_species;     /**< Number of neutral species                         */
    int* anum;         /**< Neutral species mass number                       */
    int* znum;         /**< Neutral species charge number                     */
    int* maxwellian;   /**< Is species distribution Maxwellian or
                            monoenergetic                                     */
    linint3D_data* n0; /**< Density interpolation struct for each species     */
    linint3D_data* t0; /**< Temperature intepolation struct for each species  */
} N0_3D_data;

int N0_3D_init(N0_3D_data* data,
               int n_r, real r_min, real r_max,
               int n_phi, real phi_min, real phi_max,
               int n_z, real z_min, real z_max,
               int n_species, int* anum, int* znum, int* maxwellian,
               real* density, real* temperature);
void N0_3D_free(N0_3D_data* data);
void N0_3D_offload(N0_3D_data* data);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_3D_eval_n0(real* n0, real r, real phi, real z, N0_3D_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_3D_eval_t0(real* t0, real r, real phi, real z, N0_3D_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
int N0_3D_get_n_species(N0_3D_data* ndata);
#endif
