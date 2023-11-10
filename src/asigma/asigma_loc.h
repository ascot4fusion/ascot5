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
 * @brief Local-files atomic reaction offload data
 */
typedef struct {
    int N_reac;               /**< Number of reactions                 */
    int z_1[MAX_ATOMIC];      /**< Atomic number of test particle      */
    int a_1[MAX_ATOMIC];      /**< Mass number of test particle        */
    int z_2[MAX_ATOMIC];      /**< Atomic number of bulk particle      */
    int a_2[MAX_ATOMIC];      /**< Mass number of bulk particle        */
    int N_E[MAX_ATOMIC];      /**< Size of the energy abscissa         */
    int N_n[MAX_ATOMIC];      /**< Size of the density abscissa        */
    int N_T[MAX_ATOMIC];      /**< Size of the temperature abscissa    */
    int reac_type[MAX_ATOMIC];/**< Reaction type                       */
    int offload_array_length; /**< Number of elements in offload_array */
} asigma_loc_offload_data;

/**
 * @brief Local-files atomic reaction simulation data
 */
typedef struct {
    int N_reac;                      /**< Number of reactions                 */
    int z_1[MAX_ATOMIC];             /**< Atomic number of test particle      */
    int a_1[MAX_ATOMIC];             /**< Mass number of test particle        */
    int z_2[MAX_ATOMIC];             /**< Atomic number of bulk particle      */
    int a_2[MAX_ATOMIC];             /**< Mass number of bulk particle        */
    int reac_type[MAX_ATOMIC];       /**< Reaction type                       */
    interp1D_data sigma[MAX_ATOMIC]; /**< Spline of cross-sections            */
    interp2D_data sigmav[MAX_ATOMIC];/**< Spline of rate coefficients         */
    interp3D_data BMSsigmav[MAX_ATOMIC];/**< Spline of BMS rate coefficients  */
} asigma_loc_data;

int asigma_loc_init_offload(asigma_loc_offload_data* offload_data,
                            real** offload_array);
void asigma_loc_free_offload(asigma_loc_offload_data* offload_data,
                             real** offload_array);

#pragma omp declare target
void asigma_loc_init(
    asigma_loc_data* asigma_data,
    asigma_loc_offload_data* offload_data, real* offload_array);
#pragma omp declare simd uniform(asigma_data, reac_type, z_2, a_2,\
    extrapolate)
a5err asigma_loc_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, real E_coll_per_amu,
    int reac_type, int extrapolate, asigma_loc_data* asigma_data);
#pragma omp declare simd uniform(asigma_data, reac_type, z_2, a_2,\
    extrapolate)
a5err asigma_loc_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    real E, real T_e, real T_0, real n_i, int reac_type, int extrapolate,
    asigma_loc_data* asigma_data);
#pragma omp end declare target

#endif
