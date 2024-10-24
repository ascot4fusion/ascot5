/**
 * @file N0_1D.h
 * @brief Header file for N0_1D.c
 */
#ifndef N0_1D_H
#define N0_1D_H
#include "../ascot5.h"
#include "../offload.h"
#include "../linint/linint.h" /* for 1D interpolation routines */

/**
 * @brief 1D neutral parameters on the target
 */
typedef struct {
    int n_species;     /**< Number of neutral species                         */
    int* anum;         /**< Neutral species mass number                       */
    int* znum;         /**< Neutral species charge number                     */
    int* maxwellian;   /**< Is species distribution Maxwellian or
                            monoenergetic                                     */
    linint1D_data* n0; /**< Density interpolation struct for each species     */
    linint1D_data* t0; /**< Temperature intepolation struct for each species  */
} N0_1D_data;

int N0_1D_init(N0_1D_data* data, int n_rho, real rho_min, real rho_max,
               int n_species, int* anum, int* znum, int* maxwellian,
               real* density, real* temperature);
void N0_1D_free(N0_1D_data* data);
void N0_1D_offload(N0_1D_data* data);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_1D_eval_n0(real* n0, real rho, N0_1D_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err N0_1D_eval_t0(real* t0, real rho, N0_1D_data* ndata);
DECLARE_TARGET_SIMD_UNIFORM(ndata)
int N0_1D_get_n_species(N0_1D_data* ndata);
#endif
