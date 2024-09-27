/**
 * @file asigma_loc.h
 * @brief Header file for asigma_loc.c
 */
#ifndef ASIGMALOCAL_H
#define ASIGMALOCAL_H

#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief Local-files atomic reaction simulation data
 */
typedef struct {
    int N_reac;               /**< Number of reactions             */
    int* z_1;                 /**< Atomic number of test particle  */
    int* a_1;                 /**< Mass number of test particle    */
    int* z_2;                 /**< Atomic number of bulk particle  */
    int* a_2;                 /**< Mass number of bulk particle    */
    int* reac_type;           /**< Reaction type                   */
    interp1D_data* sigma;     /**< Spline of cross-sections        */
    interp2D_data* sigmav;    /**< Spline of rate coefficients     */
    interp3D_data* BMSsigmav; /**< Spline of BMS rate coefficients */
} asigma_loc_data;


int asigma_loc_init(asigma_loc_data* data, int nreac,
                    int* z1, int* a1, int* z2, int* a2, int* reactype,
                    int* ne, real* emin, real* emax,
                    int* nn, real* nmin, real* nmax,
                    int* nT, real* Tmin, real* Tmax, real* sigma);
void asigma_loc_free(asigma_loc_data* data);
#pragma omp declare target
DECLARE_TARGET_SIMD_UNIFORM(asigma_data, reac_type, z_2, a_2,\
    extrapolate)
a5err asigma_loc_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, real E_coll_per_amu,
    int reac_type, int extrapolate, asigma_loc_data* asigma_data);
DECLARE_TARGET_SIMD_UNIFORM(asigma_data, reac_type, z_2, a_2,\
    extrapolate)
a5err asigma_loc_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    real E, real T_e, real T_0, real n_i, int reac_type, int extrapolate,
    asigma_loc_data* asigma_data);
#pragma omp declare simd uniform(asigmadata, znum, anum, nspec, extrapolate)
a5err asigma_loc_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0, int extrapolate,
    asigma_loc_data* asigmadata);
#pragma omp declare simd uniform(asigma_data, znum, anum, nion, extrapolate)
a5err asigma_loc_eval_bms(
    real* sigmav, int z_1, int a_1, real E, real mass, int nion,
    const int* znum, const int* anum, real T_e, real* n_i, int extrapolate,
    asigma_loc_data* asigma_data);
#pragma omp end declare target

#endif
