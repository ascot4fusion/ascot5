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
    int data_format;          /**< NOT IN USE                                 */
    int N_reac;               /**< Number of reactions                        */
    int offload_array_length; /**< Number of elements in offload_array, which
                                   will contain the following data for each
                                   reaction: z_1, a_1, z_2, a_2, reac_type,
                                   the number of grid points and min and max
                                   values for up to three abscissae, and an
                                   array of the cross-sections (or similar) */
} asigma_loc_offload_data;

/**
 * @brief Local-files atomic reaction simulation data
 */
typedef struct {
    int N_reac;                      /**< Number of reactions, i.e., length of
                                          below arrays                        */
    int* z_1;                        /**< Atomic number of test particle      */
    int* a_1;                        /**< Mass number of test particle        */
    int* z_2;                        /**< Atomic number of bulk particle      */
    int* a_2;                        /**< Mass number of bulk particle        */
    int* reac_type;                  /**< Reaction type                       */
    int* reac_avail;                 /**< Flag for reaction availability.
                                     NOTE: Currently without meaningful usage */
    interp1D_data* sigma;            /**< Spline of cross-sections            */
    interp2D_data* sigmav;           /**< Spline of rate coefficients         */
    interp3D_data* BMSsigmav;        /**< Spline of BMS rate coefficients     */
} asigma_loc_data;

int asigma_loc_init_offload(asigma_loc_offload_data* offload_data,
                            real** offload_array);
void asigma_loc_free_offload(asigma_loc_offload_data* offload_data,
                             real** offload_array);

#pragma omp declare target
void asigma_loc_init(asigma_loc_data* asgm_loc_data,
                     asigma_loc_offload_data* offload_data,
                     real* offload_array);
#pragma omp declare simd uniform(asgm_loc_data)
a5err asigma_loc_eval_sigma(real* sigma,
                            int z_1, int a_1,
                            int z_2, int a_2,
                            int reac_type,
                            asigma_loc_data* asgm_loc_data,
                            real E_coll_per_amu,
                            int* enable_atomic);
#pragma omp declare simd uniform(asgm_loc_data)
a5err asigma_loc_eval_sigmav(real* sigmav,
                             int z_1, int a_1, real m_1,
                             int z_2, int a_2,
                             int reac_type,
                             asigma_loc_data* asgm_loc_data,
                             real E,
                             real T_e, real T_0, real n_i,
                             int* enable_atomic);
#pragma omp end declare target

#endif
